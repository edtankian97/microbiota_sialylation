#!/bin/bash

# Define paths
FASTQ_PATH="./host_adapt_trimm_reads/"
OUTPUT_PATH="outputs_flash"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_PATH"


# Loop through all *_R1.fastq.gz files in the FASTQ_PATH
for FILE1 in "$FASTQ_PATH"/*.clean.pair.rmhost.1.fq.gz; do
  # Get the base name of the sample (remove the _R1.fastq.gz suffix)
  BASENAME=$(basename "$FILE1" .clean.pair.rmhost.1.fq.gz)

  # Define the corresponding file for the _R2.fastq pair
  FILE2="${FASTQ_PATH}/${BASENAME}.clean.pair.rmhost.1.fq.gz"

  # Check if the _R2.fastq file exists
  if [ -f "$FILE2" ]; then
    # Merging with FLASH
    flash -z -M 150 --to-stdout "$FILE1" "$FILE2" > "${OUTPUT_PATH}/${BASENAME}.merged.fastq.gz"
    # -m 20: based on the minimum merge of fastp, which is 30, and the default of FLASH, which is 10.
    # -M 150: avoids the WARNING from FLASH saying that many sequences overlap more than the parameter, but it didn't change the merge percentage
  else
    echo "Pair for file $FILE1 not found: $FILE2" >&2
  fi
done
