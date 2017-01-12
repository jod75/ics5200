#!/bin/bash

hdfs dfsadmin -safemode enter

$SPARK_HOME/sbin/stop-all.sh

$HADOOP_HOME/sbin/stop-dfs.sh
$HADOOP_HOME/sbin/stop-yarn.sh
