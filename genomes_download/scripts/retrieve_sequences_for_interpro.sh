#! /bin/bash


for file in ../proteins/proteins_sialylation/*for_Interpro.tsv; do

BASE_FILE=$(basename "$file" .tsv)
seqkit grep -f "$file" -n *.faa > "${BASE_FILE}_retrieved_now"
done
