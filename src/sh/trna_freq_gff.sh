#!/bin/bash

for f in /data3/clinical_isolates_sq/prokka_gff/*.gff
do
  FNAME=`basename ${f%%.*}`;
  awk -F "\t" '$3 == "tRNA" {print $9}' $f |  awk -F ";" '{split($5, words, "="); print words[2]}' | sort | uniq -c | sort -nr | nl > /home/dst20/dev/ribotracker/dat/tmp/${FNAME}_tRNA_gene_count.txt;
  echo "Processing strain ${FNAME}";
  echo "Writing into /home/dst20/dev/ribotracker/dat/tmp/${FNAME}_tRNA_gene_count.txt";
done
