#!/usr/bin/env bash
PYSPARK_DRIVER_PYTHON="jupyter" PYSPARK_DRIVER_PYTHON_OPTS="notebook" pyspark --master yarn --verbose --packages graphframes:graphframes:0.3.0-spark2.0-s_2.11
