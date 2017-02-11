#!/bin/bash

USERNAME=hduser
HOSTS="hadoop1 hadoop2 hadoop3 hadoopa hadoopb"
SCRIPT="sudo apt-get install ncbi-blast+; mkdir /home/hduser/proteins; rm -r /home/hduser/proteins/*"
for HOSTNAME in ${HOSTS} ; do
    echo ${HOSTNAME}
    ssh -t -o StrictHostKeyChecking=no -l ${USERNAME} ${HOSTNAME} "${SCRIPT}"
    scp /home/hduser/Lab/prosim/proteins.* ${USERNAME}@${HOSTNAME}:/home/hduser/proteins/
    scp /home/hduser/Lab/702blast.sh ${USERNAME}@${HOSTNAME}:/home/hduser/proteins/702blast.sh
done

