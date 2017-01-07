#!/bin/bash

cd ~
$SPARK_HOME/sbin/stop-all.sh

$HADOOP_HOME/sbin/stop-yarn.sh
$HADOOP_HOME/sbin/stop-dfs.sh

USERNAME=hduser
HOSTS="hadoop2 hadoop3"
SCRIPT="sudo shutdown -P now"
for HOSTNAME in ${HOSTS} ; do
    echo ${HOSTNAME}
    ssh -t -o StrictHostKeyChecking=no -l ${USERNAME} ${HOSTNAME} "${SCRIPT}"
done

ssh -t -o StrictHostKeyChecking=no joseph@192.168.151.11 "sudo shutdown -P now"

echo Goodbye... powering down myself...
sudo shutdown -P now
