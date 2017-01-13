#!/usr/bin/env python

##################################################################################
# engine.py
from pyspark.sql.types import *
from pyspark.sql.functions import col, desc, concat, lit
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
    def __setup(self):
        """ Performs start-up routines:
            1. Loads data from ChEMBL MySQL and store it in a TSV file on Hadoop
            2. Take samples for testing:
                a. Ligands
                b. Proteins
            3. Create DataFrames:
                a. Ligands:
                    - Ligands = unique list of Ligands
                    - LigandsBindingsTest = test ligands known bindings used for evaluation
                    - LigandsBindingsSample = sample ligands known bindings
                b. Proteins:
                    - Proteins = unique list of Proteins
                    - ProteinsBindingsTest = test proteins known bindings used for evaluation
                    - ProteinsBindingsSample = sample proteins known bindings
            4. Create local blast database
        """

        # Step 1 - Get raw data
        chemblhelper = ChEMBLHelper()

        hdfsServer = "http://hadoop1:50070" # hdfs path
        localHome = "/home/hduser/Lab"
        hdfsHome = "/user/hduser/ics5200"
        datasetCount = 100000 # dataset count of bindings from ChEMBL

        self.sparkFastaFile = "/home/hduser/Lab/chembl" + \
            str(datasetCount) + \
            ".fasta" # file created by spark task
        self.localFasta = "/home/hduser/Lab/proteinbank/chembl" + str(datasetCount) + ".fasta"
        self.blastOutFile = "/home/hduser/Lab/blastResult.xml"

        self.databankFilename = chemblhelper.createDataBank(hdfsServer,
                                                            localHome,
                                                            hdfsHome,
                                                            datasetCount)

        # Step 2 and 3 - Sampling
        self.__sampleData(0.01, 1) # 1% used for testing

        # Step 4 - Create Blast DB
        self.__createLocalBlastDb()

    def __sampleData(self, testPercentage=0.01, seeding=None):
        """  Splits raw data into test and known sets
        """
        # dataRDD will contain all data imported from MySQL
        dataRDD = self.sc.textFile(self.databankFilename).map(lambda line: line.split('\t'))

        # convert each line (currently an list) to a tuple.  This makes it easier to manipulate
        # the data, especially to convert to DataFrames
        dataRDD = dataRDD.map(lambda l: tuple(l))

        ### unique lists
        # Ligands - molRergNo
        testLigandsList = dataRDD.map(lambda t: long(t[2])).distinct().sample(withReplacement=False, seed=seeding, fraction=testPercentage)
        # Proteins - compId (which is long rather than accession, which is a string)
        testProteinsList = dataRDD.map(lambda t: long(t[8])).distinct().sample(withReplacement=False, seed=seeding, fraction=testPercentage)

        # broadcast value as a list of lonf numbers so that mapping is faster
        testLigandsListBV = self.sc.broadcast(testLigandsList.collect())
        testProteinsListBV = self.sc.broadcast(testProteinsList.collect())

        ### samples
        # Test sets
        ligandsBindingsTestRDD = dataRDD.filter(lambda t: long(t[2]) in testLigandsListBV.value)
        proteinsBindingsTestRDD = dataRDD.filter(lambda t: long(t[8]) in testProteinsListBV.value)

        # Known sets
        ligandsBindingsKnownRDD = dataRDD.subtract(ligandsBindingsTestRDD)
        proteinsBindingsKnownRDD = dataRDD.subtract(proteinsBindingsTestRDD)

        ### DataFrames
        # Schemas
        ligandsSchema = StructType([StructField("molregno", IntegerType(), False),
                                    StructField("canonical_smiles", StringType(), False),
                                    StructField("mol_pref_name", StringType(), True)])

        proteinsSchema = StructType([StructField("component_id", LongType(), False),
                                     StructField("accession", StringType(), True),
                                     StructField("sequence", StringType(), False),
                                     StructField("prot_pref_name", StringType(), True),
                                     StructField("prot_short_name", StringType(), True)])

        bindingsSchema = StructType([StructField("row_id", LongType(), False),
                                     StructField("assay_id", LongType(), False),
                                     StructField("molregno", LongType(), False),
                                     StructField("std_relation", StringType(), True),
                                     StructField("std_value", DecimalType(), True),
                                     StructField("std_units", StringType(), True),
                                     StructField("std_type", StringType(), True),
                                     StructField("pchembl_value", DecimalType(), True),
                                     StructField("component_id", LongType(), False)])

        # Unique Lists
        ligandsRDD = dataRDD.map(lambda t: (long(t[2]), str(t[11]), str(t[12]))).distinct()
        proteinsRDD = dataRDD.map(lambda t: (long(t[8]), str(t[9]), str(t[10]), str(t[13]), str(t[14]))) \
                             .distinct()
        self.ligands = self.sqlContext.createDataFrame(ligandsRDD, ligandsSchema)
        self.proteins = self.sqlContext.createDataFrame(proteinsRDD, proteinsSchema)

        # Test data sets
        self.ligandsBindingsTest = self.sqlContext.createDataFrame(ligandsBindingsTestRDD.map(
            lambda t: (long(t[0].rstrip(".0")),
                       long(t[1]),
                       long(t[2]),
                       str(t[3]),
                       PythonHelper.getDecimal(t[4]),
                       str(t[5]),
                       str(t[6]),
                       PythonHelper.getDecimal(t[7]),
                       long(t[8]))), bindingsSchema)

        self.proteinsBindingsTest = self.sqlContext.createDataFrame(proteinsBindingsTestRDD.map(
            lambda t: (long(t[0].rstrip(".0")),
                       long(t[1]),
                       long(t[2]),
                       str(t[3]),
                       PythonHelper.getDecimal(t[4]),
                       str(t[5]),
                       str(t[6]),
                       PythonHelper.getDecimal(t[7]),
                       long(t[8]))), bindingsSchema)

        # Known data sets
        self.ligandsBindingsKnown = self.sqlContext.createDataFrame(ligandsBindingsKnownRDD.map(
            lambda t: (long(t[0].rstrip(".0")),
                       long(t[1]),
                       long(t[2]),
                       str(t[3]),
                       PythonHelper.getDecimal(t[4]),
                       str(t[5]),
                       str(t[6]),
                       PythonHelper.getDecimal(t[7]),
                       long(t[8]))), bindingsSchema)

        self.proteinsBindingsKnown = self.sqlContext.createDataFrame(proteinsBindingsKnownRDD.map(
            lambda t: (long(t[0].rstrip(".0")),
                       long(t[1]),
                       long(t[2]),
                       str(t[3]),
                       PythonHelper.getDecimal(t[4]),
                       str(t[5]),
                       str(t[6]),
                       PythonHelper.getDecimal(t[7]),
                       long(t[8]))), bindingsSchema)

    def __createLocalBlastDb(self):
        """ Create BLAST local db
        """

        # clean up before saving file to disk
        shutil.rmtree(self.sparkFastaFile, ignore_errors=True)

        # manipulate raw data rdd and create FASTA file
        blastDb = self.proteinsBindingsKnown \
                    .join(self.proteins,
                          self.proteinsBindingsKnown.component_id == self.proteins.component_id) \
                    .select("accession", "prot_pref_name", "sequence") \
                    .distinct() \
                    .rdd \
                    .map(lambda t: str(">ebl|" + t[0] + "| " + t[1] + "\r" + "\r".join(re.findall(".{1,80}",t[2]))))
        
        dump = blastDb.collect()
        myfile = open(self.localFasta, "w")
        for item in dump:
            myfile.write("%s\r" % item)            
        
        # create local blast db
        process = subprocess.Popen(shlex.split(
            "makeblastdb -in %s -parse_seqids -dbtype prot" % os.path.basename(self.localFasta)),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,  
            cwd=os.path.dirname(self.localFasta),
            shell=False)

        out, err = process.communicate()

        PythonHelper.writeToJupyterConsole(">Engine init makeblastdb out: " + out)
        PythonHelper.writeToJupyterConsole(">Engine init makeblastdb err: " + err)

    def getTestLigandsDS(self):
        return self.ligandsBindingsTest.collect()

    def getTestProteinsDS(self):
        return self.proteinsBindingsTest.collect()

    def getSmiles(self, molRegNo):
        """ Get SMILES representation of a molecule.

            Args:
                molRegNo: long, molRegNo in ChEMBL                

            Returns:
                string: SMILES representation.
        """

        # there should be only one entry and one field
        return self.ligands.filter(col("molregno") == long(molRegNo)).select("canonical_smiles").collect()[0][0]

    def getLigandsTestTargets(self, molRegNo):
        """ Returns the bindings in test data set.

            Args:
                molRegNo: long
            
            Returns:
                List[long]: a list of component_id
        """

        return self.ligandsBindingsTest.filter(col("molregno") == molRegNo).select("component_id").distinct().orderBy("component_id").collect()

    def doLigandExperiment(self, molRegNo, molHelper = MoleculeHelper, similarityThreshold = 0.5):
        """ Runs a ligand experiment

            Args:
                molRegNo: long

            Returns:
                List of known bindings
        """
                
        querySmiles = self.getSmiles(molRegNo)
        queryLigand = dict()
        queryLigand.update({0: querySmiles})
        queryRDD = self.sc.parallelize(queryLigand).map(lambda k:(k, molHelper(queryLigand[k])))
        ligandsRDD = self.ligandsBindingsKnown.join(self.ligands, self.ligandsBindingsKnown.molregno == self.ligands.molregno) \
                                              .select(self.ligands.molregno, self.ligands.canonical_smiles) \
                                              .rdd.map(lambda (k, v):(k, molHelper(v)))                                              
        simRDD = ligandsRDD.cartesian(queryRDD) \
                           .map(lambda ((k1,v1),(k2,v2)): (k1, float(v1.similarity(v2)))) \
                           .filter(lambda (k1, v): v >= similarityThreshold)

        simSchema = StructType([StructField("molregno", LongType(), False),                                
                                StructField("similarity", FloatType(), False)])
        sim = self.sqlContext.createDataFrame(simRDD, simSchema)
        
        return self.ligandsBindingsKnown.join(sim, self.ligandsBindingsKnown.molregno == sim.molregno).orderBy(desc("similarity")).collect()
    
    def __init__(self, sc, dataset_path):
        """Init the recommendation engine given a Spark context and a dataset path
        """

        logger.info("Starting up the ICS5200 Engine: ")

        self.sc = sc
        self.sqlContext = SQLContext(sc)

        # Load ratings data for later use
        logger.info("Loading data...")

        self.__setup()

