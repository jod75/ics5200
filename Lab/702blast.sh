#!/bin/sh
while read LINE; do
  q=${LINE}
  read LINE
  s=${LINE}
  m=$q'\n'$s
  echo $m | blastp -db /home/hduser/proteins/proteins.fasta -outfmt 10 -query - 2>&1
done
