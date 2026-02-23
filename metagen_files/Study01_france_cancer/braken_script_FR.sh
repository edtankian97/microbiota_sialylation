#! /bin/bash

OUT=bracken_out

mkdir -p $OUT

for file in kraken_out/*.report; do

BASE=$(basename $file .report)

bracken -d kraken2_DB/ -i $file -l S -o ${OUT}/${BASE}_bracken_out.tsv;

done
