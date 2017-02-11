#!/bin/sh

hdfs dfsadmin -safemode leave
hdfs fsck /
