#!/usr/bin/env python

##################################################################################
# engine.py
from pyspark.sql.types import *
from pyspark.sql.functions import col, desc
from pyspark.sql import SQLContext

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import MACCSkeys

from moleculehelper import *
from chemblhelper import ChEMBLHelper
from pythonhelper import *
from hdfshelper import HDFSHelper

import os.path
import logging
import json

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class ICS5200Engine(object):
    def __seteup(self):
        # get small dataset from ChEMBL
        chemblhelper = ChEMBLHelper()

        hdfsServer = "http://hadoop1:50070"                          # hdfs path
        datasetCount = 100000                                        # dataset count of bindings from ChEMBL
        datasetTSVFilename = "sample" + str(datasetCount) + ".tsv"   # file with sample data
        hdfsDatasetFilename = os.path.join("/user/hduser", datasetTSVFilename)

        # check if the dataset file exists in hdfs and if it does not, then load data from ChEMBL database
        if not HDFSHelper.fileExists(hdfsServer, hdfsDatasetFilename):
            # get data from ChEMBL and stores the file to local dfs
            chemblhelper.saveBindingsToTSV(datasetTSVFilename, datasetCount) 
            # upload data to hdfs so that it is accessbile from all cluster worker nodes
            HDFSHelper.putFile(hdfsServer, datasetTSVFilename, hdfsDatasetFilename)     
                    
        dataRDD = self.sc.textFile(hdfsDatasetFilename).map(lambda line: line.split('\t'))
        
        # convert each line (currently an list) to a tuple.  This makes it easier to manipulate 
        # the data, especially to convert to DataFrames
        dataRDD = dataRDD.map(lambda l: tuple(l))
        
        moleculesSample = dataRDD.map(lambda t: long(t[2])).distinct().sample(withReplacement=False, seed=1, fraction=0.01)
        
        # broadcast value as a list of lonf numbers so that mapping is faster
        moleculesSampleBV = self.sc.broadcast(moleculesSample.collect())
        
        # test data is reserved as unseen knowledge
        self.testData = dataRDD.filter(lambda t: long(t[2]) in moleculesSampleBV.value)
        
        # sample data is the known dataset
        sampleData = dataRDD.subtract(self.testData)
        
        # each line need to be converted to a tuple so that later it is converted into a DF
        moleculesRDD = sampleData.map(lambda t: (long(t[2]),str(t[11]))).distinct()

        # binding schema
        moleculesSchema = StructType([StructField("molregno", IntegerType(), False),
                                      StructField("canonical_smiles", StringType(), False)])

        # convert RDD to DataFrame - faster and more memory efficient
        self.molecules = self.sqlContext.createDataFrame(moleculesRDD, moleculesSchema)

        bindingsRDD = sampleData.map(lambda t: (long(t[0].rstrip(".0")),
                                        long(t[1]),                                      
                                        long(t[2]),
                                        str(t[3]),
                                        PythonHelper.getDecimal(t[4]),
                                        str(t[5]),
                                        str(t[6]),
                                        PythonHelper.getDecimal(t[7]),
                                        long(t[8]))
                            ) # each line need to be converted to a tuple so that later it is converted into a DF

        # binding schema
        bindingsSchema = StructType([StructField("row_id", LongType(), False),
                                     StructField("assay_id", LongType(), False),                            
                                     StructField("molregno", LongType(), False),
                                     StructField("std_relation", StringType(), True),
                                     StructField("std_value", DecimalType(), True),
                                     StructField("std_units", StringType(), True),
                                     StructField("std_type", StringType(), True),
                                     StructField("pchembl_value", DecimalType(), True),
                                     StructField("component_id", LongType(), False)])

        # convert RDD to DataFrame - faster and more memory efficient
        self.bindings = self.sqlContext.createDataFrame(bindingsRDD, bindingsSchema)
        
        # each line need to be converted to a tuple so that later it is converted into a DF
        proteinsRDD = sampleData.map(lambda t: (long(t[8]),str(t[9]),str(t[10]))).distinct()

        # binding schema
        proteinsSchema = StructType([StructField("component_id", LongType(), False),
                                     StructField("accession", StringType(), True),
                                     StructField("sequence", StringType(), False)])

        # convert RDD to DataFrame - faster and more memory efficient
        self.proteins = self.sqlContext.createDataFrame(proteinsRDD, proteinsSchema)

    # this must run on the main thread
    def __findSimilarMolecules(self, querySmiles, knownMolecules, molHelper = MoleculeHelper, similarityThreshold = 0.85):
        """ Returns an RDD with similar molecules.
        """
        
        # step 1 - create a molecule helper class for each molecule, this will take
        #          more memory but will increase computation efficiency
        queryMol = dict()
        queryMol.update({0: querySmiles})
        queryRDD = self.sc.parallelize(queryMol).map(lambda k:(k, molHelper(queryMol[k])))
        mols = knownMolecules.rdd.map(lambda (k, v):(k, molHelper(v)))

        sm = mols.cartesian(queryRDD) \
                 .map(lambda ((k1,v1),(k2,v2)): (k1, float(v1.similarity(v2)))) \
                 .filter(lambda (k1, v): v >= similarityThreshold)    
                
        simSchema = StructType([StructField("molregno", LongType(), False),
                                #StructField("queryMol", LongType(), False),
                                StructField("similarity", FloatType(), False)])

        sim = self.sqlContext.createDataFrame(sm, simSchema)
        
        return sim
        
    def __getBindings(self, similarMolecules):
        return self.bindings.join(similarMolecules, self.bindings.molregno == similarMolecules.molregno)
        
    def __experiment(self, experimentsData, moleculeIndex, knownMolecules, molHelper = MoleculeHelper, similarityThreshold = 0.85):
        selectedMolecules = experimentsData.map(lambda t: (long(t[2]), t[11])).distinct().sortBy(lambda x: x[0])
        selectedMoleculeMolRegNo = selectedMolecules.collect()[moleculeIndex][0]
        selectedMoleculeSMILES = selectedMolecules.collect()[moleculeIndex][1]
        
        logger.info( "MolRegNo = " + str(selectedMoleculeMolRegNo))
        logger.info( "SMILES = " + selectedMoleculeSMILES)
        
        expData = experimentsData.filter(lambda t: long(t[2]) == selectedMoleculeMolRegNo).map(lambda m: m[8]).distinct().collect()
        
        logger.info( "Target Componet IDs:")
        for compId in expData:
            logger.info( "  " + str(compId))
        
        sm = self.__findSimilarMolecules(selectedMoleculeSMILES, knownMolecules, molHelper, similarityThreshold)
        return self.__getBindings(sm)  
        
    def getMoleculeInfo(self, userId, moleculeId):        
        test = self.__experiment(self.testData, moleculeId, self.molecules, similarityThreshold = 0.5).orderBy(desc("similarity")).collect()
        return test
        
    def getSmiles(self, molregno):
        return self.molecules.filter(col("molregno") == molregno).take(1)
        
    def getTestSmiles(self, molIndex):
        return self.testData.map(lambda t: (t[2], t[11])).distinct().collect()[molIndex][1]
        
    def getTestMolRegNos(self):
        x = self.testData.map(lambda t: t[2]).distinct().count()
        return range(0, x)

    def __init__(self, sc, dataset_path):
        """Init the recommendation engine given a Spark context and a dataset path
        """

        logger.info("Starting up the ICS5200 Engine: ")

        self.sc = sc      
        self.sqlContext = SQLContext(sc)

        # Load ratings data for later use
        logger.info("Loading data...")
        
        self.__seteup()

