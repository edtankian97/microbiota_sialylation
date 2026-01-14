#! /bin/bash

OUT_PATH="./output_data/Interpro_results/"

mkdir -p $OUT_PATH

for file in ./output_data/Interpro_analysis/*retrieved; do

BASENAME=$(basename $file)

../../genomes_download/interproscan-5.76-107.0/interproscan.sh -i $file -b $OUT_PATH/${BASENAME}_results -f tsv -appl CDD,COILS,Gene3D,HAMAP,MobiDBLite,PANTHER,Pfam,PIRSF,PRINTS,SFLD,SMART,SUPERFAMILY,TIGRFAM;
done
