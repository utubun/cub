#!/bin/bash

PATH = ./dat/raw/fna

for FILES in /data3/clinical_isolates_sq/DNA-seq/final_data_all_clinicals/prokka_annotation/fna_files/CH25*.fna
do
  FNAME=`basename ${f%%.*}`;
  tRNAscan-SE -B -L -H -D -o /home/dst20/dev/ribotracker/dat/trnascan/${FNAME}.out -f /home/dst20/dev/ribotracker/dat/trnascan/${FNAME}.struct  -a /home/dst20/dev/ribotracker/dat/tmp/${FNAME}.fa -m /home/dst20/dev/ribotracker/dat/tmp/${FNAME}.stats -l /home/dst20/dev/ribotracker/dat/tmp/${FNAME}.log -d -Q $f;
  echo "Processing strain ${FNAME}";
  echo "Writing into /home/dst20/dev/ribotracker/dat/trnascan/ under ${FNAME} identifier";
done
