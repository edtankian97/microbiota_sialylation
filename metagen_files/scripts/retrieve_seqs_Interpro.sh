#!/bin/bash


for file in colon*; do

BASE_FILE=$(basename "$file" .tsv)
seqkit grep -f "$file" -n *.faa > "${file}_retrieved_now"
done
