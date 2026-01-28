#! /bin/bash

OUT_PATH="../Interpro_analysis/Interpro_results/"

mkdir -p $OUT_PATH

for file in ../Interpro_analysis/*_retrieved_final; do

BASENAME=$(basename $file)

../interproscan-5.76-107.0/interproscan.sh -i $file -b $OUT_PATH/${BASENAME}_results -f tsv -appl CDD,COILS,Gene3D,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,SFLD,SMART,SUPERFAMILY,TIGRFAM;
done
