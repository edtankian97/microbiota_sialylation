#! /bin/bash


for file in ../proteins/proteins_sialylation/*for_Interpro.tsv; do

BASE_FILE=$(basename "$file" .tsv)
seqkit grep -f "$file" -n *.faa > "${file}_retrieved_now"
done
