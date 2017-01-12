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

from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML

from moleculehelper import *
from chemblhelper import ChEMBLHelper
from pythonhelper import *
from hdfshelper import HDFSHelper

import os.path
import logging
import json
import subprocess
import shutil
import shlex
import re

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
        sparkFastaFile = "/home/hduser/Lab/chembl"+ str(datasetCount) +".fasta" # file created by spark task
        self.localFasta = "/home/hduser/Lab/proteinbank/chembl" + str(datasetCount) + ".fasta"
        self.blastOutFile = "/home/hduser/Lab/blastResult.xml"

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
        
        # this RDD will is used to get unique molecules from test data
        self.testMoleculeData = self.testData.map(lambda t: (long(t[2]),str(t[11]))).distinct()
        
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
        
        # create BLAST local db
        # clean up before saving file to disk
        if not os.path.exists(self.localFasta):            
            shutil.rmtree(sparkFastaFile, ignore_errors=True)
            # manipulate raw data rdd and create FASTA file
            sampleData.map(lambda line: (line[9], line[10])) \
                        .distinct() \
                        .map(lambda t: str(">ebl|" + t[0] + "|\r" + "\r".join(re.findall(".{1,80}",t[1])))) \
                        .coalesce(1) \
                        .saveAsTextFile("file://" + sparkFastaFile)
            os.rename(os.path.join(sparkFastaFile, "part-00000"), self.localFasta)
            
            # create local blast db
            process = subprocess.Popen(shlex.split("makeblastdb -in %s -parse_seqids -dbtype prot" % os.path.basename(self.localFasta)),                           
                                       stdout = subprocess.PIPE,
                                       stderr = subprocess.PIPE,  
                                       cwd = os.path.dirname(self.localFasta),
                                       shell = False)

            out, err = process.communicate()

            PythonHelper.writeToJupyterConsole(">Engine init makeblastdb out: " + out)
            PythonHelper.writeToJupyterConsole(">Engine init makeblastdb err: " + err)
        

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
        
    def __experiment(self, experimentsData, molRegNo, knownMolecules, molHelper = MoleculeHelper, similarityThreshold = 0.85):
        #selectedMoleculeMolRegNo = self.testMoleculeData.collect()[moleculeIndex][0]
        selectedMoleculeSMILES = self.getTestSmiles(molRegNo)
        
        logger.info( "MolRegNo = " + str(molRegNo))
        logger.info( "SMILES = " + selectedMoleculeSMILES)
        
        expData = experimentsData.filter(lambda t: long(t[2]) == molRegNo).map(lambda m: m[8]).distinct().collect()
        
        logger.info( "Target Componet IDs:")
        for compId in expData:
            logger.info( "  " + str(compId))
        
        sm = self.__findSimilarMolecules(selectedMoleculeSMILES, knownMolecules, molHelper, similarityThreshold)
        return self.__getBindings(sm)  
     
    # run experiment
    def doExperiment(self, molRegNo):
        """ Search for molecules similar to molRegNo.
            
            Args:
                molRegNo: The molRegNo in ChEMBL database for a molecule in test data set.
                
            Returns:
                List: An RDD with similar molecules and their bindings.
        """      
        
        #test = self.__experiment(self.testData, molRegNo, self.molecules, molHelper = MoleculeMACCSHelper, similarityThreshold = 0.85).orderBy(desc("similarity")).collect()
        test = self.__experiment(self.testData, molRegNo, self.molecules, similarityThreshold = 0.5).orderBy(desc("similarity")).collect()
        return test
        
    # get SMILES for a molecule from Sample Data
    def getSmiles(self, molRegNo):
        """ Get SMILES for a molecule in Sample data set.
        
            Args:
                molRegNo: The molRegNo in ChEMBL database for a molecule in Sample data set for which SMILES string sequence is required.
                
            Returns:
                String: Canonical SMILES representation.
        """
        
        return self.molecules.filter(col("molregno") == long(molRegNo)).take(1)[0][1]
        
    # get SMILES for a molecule in Test Data
    def getTestSmiles(self, molRegNo):
        """  Get SMILES for a molecule in Test data set.
        
            Args:
                molRegNo: The molRegNo in ChEMBL database for a molecule in Test data set for which SMILES string sequence is required.
                
            Returns:
                String: Canonical SMILES representation.
        """
        
        return self.testMoleculeData.filter(lambda t: t[0] == long(molRegNo)).take(1)[0][1]
        
    # get a list of targets
    def getTestTargets(self, molRegNo):  
        """  Get list of compound target id for a molecule in Test data set.
        
            Args:
                molRegNo: The molRegNo in ChEMBL database for a molecule in Test data set for which targets list is required.
                
            Returns:
                List: Compound target ids.
        """
        
        return self.testData.filter(lambda t: t[2] == molRegNo).map(lambda m: m[8]).distinct().collect()
        
    # get a data set
    def getTestDataset(self):  
        """ Gets all data in Test data set.
        
            Returns:
                List: All data in test data set.
        """
        
        return self.testData.collect()
        
    # get the molRegNos in Test Data
    def getTestMolRegs(self):
        """ Gets list of molRegNos in test data set.
        
            Returns:
                List: molRegNos in test data.
        """
        
        return self.testMoleculeData.map(lambda t: long(t[0])).collect()
        
    def doProteinExperiment(self, compId):
        query_seq = self.testData.filter(lambda t: t[8] == compId).distinct().map(lambda m: m[10]).collect()[0]
        blastp_cline = NcbiblastxCommandline(cmd = "blastp",
                                     db = self.localFasta,
                                     evalue = 0.01,
                                     outfmt = 5,
                                     remote = False,
                                     out = self.blastOutFile)

        (out, err) = blastp_cline(stdin = query_seq)
        PythonHelper.writeToJupyterConsole(">Engine doProteinExperiment out: " + out)
        PythonHelper.writeToJupyterConsole(">Engine doProteinExperiment err: " + err)
        
        # create matching proteins list.
        # each list entry is a tuple in the format:
        #  (Accession, (Hit, High-Scoring-Pair expect value))
        #  nested tuples are needed to join rdds
        result_handle = open(self.blastOutFile)
        blast_record = NCBIXML.read(result_handle)
        matches = []
        for alignment in blast_record.alignments:   
            i = 1
            for hsp in alignment.hsps:
                matches.append((alignment.accession, (i, hsp.expect)))
                i = i + 1 
        
        sm = self.sc.parallelize(matches)       
        uniqueAccessions = sm.map(lambda t: t[0]).distinct().collect()
        # this must work on dataframe not RDD
        accessionsCompId = self.proteins.rdd.filter(lambda t: t[1] in uniqueAccessions).map(lambda p: (p[1], p[0])).distinct()
        simSchema = StructType([StructField("accession", StringType(), False),
                                StructField("hit", IntegerType(), False),
                                StructField("similarity", FloatType(), False),
                                StructField("compId", LongType(), False)])

        sim = self.sqlContext.createDataFrame(sm.join(accessionsCompId).map(lambda t: (t[0], t[1][0][0], t[1][0][1], long(t[1][1]))), simSchema)        
        return self.bindings.join(sim, self.bindings.component_id == sim.compId).collect()
        

    def __init__(self, sc, dataset_path):
        """Init the recommendation engine given a Spark context and a dataset path
        """

        logger.info("Starting up the ICS5200 Engine: ")

        self.sc = sc      
        self.sqlContext = SQLContext(sc)

        # Load ratings data for later use
        logger.info("Loading data...")
        
        self.__seteup()

