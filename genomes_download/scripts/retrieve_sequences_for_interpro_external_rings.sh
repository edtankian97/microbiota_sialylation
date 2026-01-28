#! /bin/bash

for file in ../proteins/proteins_sialylation/final_complete_sialylation/*_ID_final.tsv; do

BASE_FILE=$(basename "$file" .tsv)
seqkit grep -f "$file" -n *.faa > "${BASE_FILE}_retrieved_final"
mv "${BASE_FILE}_retrieved_final" ../Interpro_analysis
done
