#!/bin/bash 

# Dependencies: MEGAHIT
# Install:
#   conda install -c bioconda megahit

FASTQ_PATH="./clean_reads_FR/"
MERGED_PATH="./outputs_flash/"
OUTPUT_PATH="./megahit_outputs/"

mkdir -p "$OUTPUT_PATH"

for FILE1 in "$FASTQ_PATH"/ERR47*_aligned_R1.fastq.gz; do

  BASENAME=$(basename "$FILE1" _aligned_R1.fastq.gz)

  FILE2="${FASTQ_PATH}/${BASENAME}_aligned_R2.fastq.gz"

  FILE3="${MERGED_PATH}/${BASENAME}_aligned.merged.fastq"

  megahit \
    -1 "$FILE1" \
    -2 "$FILE2" \
    -r "$FILE3" \
    -o "$OUTPUT_PATH/${BASENAME}" \
    -t 12

done
