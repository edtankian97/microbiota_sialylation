#!/bin/bash

# file that contains the list of files to move
for FILE1 in all_*_colon02_output_coverage.tsv_filtered_ID; do

# create output directory based on the file name (without path)
OUT_DIR="$(basename "$FILE1" _output_coverage.tsv_filtered_ID)"

mkdir -p "$OUT_DIR"

while read -r file; do
  mv "$file" "$OUT_DIR/"
done < "$FILE1"

done
