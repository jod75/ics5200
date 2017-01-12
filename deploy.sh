#!/bin/bash

USERNAME=hduser
HOSTS="hadoop2 hadoop3"
#SCRIPT="sudo apt-get install python-rdkit librdkit1 rdkit-data"
SCRIPT="sudo apt-get install python-pip; sudo pip install biopython"
for HOSTNAME in ${HOSTS} ; do
    echo ${HOSTNAME}
    ssh -t -o StrictHostKeyChecking=no -l ${USERNAME} ${HOSTNAME} "${SCRIPT}"
done

