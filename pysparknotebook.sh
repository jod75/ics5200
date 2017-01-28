#!/usr/bin/env bash

# loads jupyter notebook service and use pyspark kernel
# loads also the graphframes package
# change --master yarn to --master local
# see http://spark.apache.org/docs/latest/submitting-applications.html 
PYSPARK_DRIVER_PYTHON="jupyter" PYSPARK_DRIVER_PYTHON_OPTS="notebook" pyspark --master yarn --verbose --packages graphframes:graphframes:0.3.0-spark2.0-s_2.11
