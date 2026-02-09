#!/bin/bash


for file in ./bins_paired/all_prokka/*for_Interpro_france.tsv; do

BASE_FILE=$(basename "$file" .tsv)
seqkit grep -f "$file" -n /bins_paired/all_prokka/*.faa > "${file}_retrieved_now"
done
