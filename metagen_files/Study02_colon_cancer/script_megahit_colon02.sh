#!/bin/bash
# Dependencies: MEGAHIT
# Install:
#   conda install -c bioconda megahit

#define pathways
FASTQ_PATH="./host_adapt_trimm_reads/"
MERGED_PATH="./outputs_flash_colon_china/"   #extra file with merged sequences for better assembling
OUTPUT_PATH="megahit_outputs_colon_china"

#Create output directory if not there
mkdir -p "$OUTPUT_PATH"

#Loop to proccess all files in FASTQ_PATH
for FILE1 in "$FASTQ_PATH"/*.clean.trim.rmhost.1.fq.gz; do
  # Obtain the basename of the sample (remove sufix _filtered_aligned_R1. fastq)
  BASENAME=$(basename "$FILE1" .clean.trim.rmhost.1.fq.gz)

  # Define the pair file _filtered_aligned_R2.fastq
  FILE2="${FASTQ_PATH}/${BASENAME}.clean.trim.rmhost.2.fq.gz"

  #Define the merged file
  FILE3="${MERGED_PATH}/${BASENAME}.merged.fastq.gz"

  #Assembling with MEGAHIT
  megahit -1 "$FILE1" -2 "$FILE2" -r "$FILE3" -o "$OUTPUT_PATH/${BASENAME}" -t 5
done
