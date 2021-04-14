#!/bin/bash

SRCPATH=$(ls ./dat/raw/aln/fna/CH5*)

for FILE in ${SRCPATH}
do
  FNAME=`basename -s .fna ${FILE}`;

  OUTPATH="./dat/out/trs/output/${FNAME}.out";
  STRPATH="./dat/out/trs/struct/${FNAME}.struct";
  FASPATH="./dat/out/trs/fasta/${FNAME}.fa";
  BEDPATH="./dat/out/trs/bed/${FNAME}.bed";
  STSPATH="./dat/out/trs/stats/${FNAME}.stats";
  LOGPATH="./dat/out/trs/log.log";

  tRNAscan-SE -B -L -H -D -o ${OUTPATH} -f ${STRPATH} -a ${FASPATH} -b ${BEDPATH} -m ${STSPATH} -l ${LOGPATH} -d -Q $FILE;
  echo "Processing strain ${FNAME}";
done
