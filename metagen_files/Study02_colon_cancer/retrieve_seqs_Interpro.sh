#!/bin/bash


for file in colon02*; do

BASE_FILE=$(basename "$file" .tsv)
seqkit grep -f "$file" -n bins_paired/all_prokka/*.faa > "${file}_retrieved_now"
done
