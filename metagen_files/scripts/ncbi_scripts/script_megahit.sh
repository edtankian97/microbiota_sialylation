#! /bin/bash 

# Dependencies: MEGAHIT
# Install:
#   conda install -c bioconda megahit

#define pathways
FASTQ_PATH="./clean_reads/"
MERGED_PATH="./outputs_flash/"   #extra file with merged sequences for better assembling
OUTPUT_PATH="./megahit_outputs/"

#Create output directory if not there
mkdir -p "$OUTPUT_PATH"


#Loop to proccess all files in FASTQ_PATH
for FILE1 in "$FASTQ_PATH"/*_; do
  # Obtain the basename of the sample (remove sufix _filtered_aligned_R1. fastq)
  BASENAME=$(basename "$FILE1" _filtered_aligned_R1.fastq.gz)

  # Define the pair file _filtered_aligned_R2.fastq
  FILE2="${FASTQ_PATH}/${BASENAME}_filtered_aligned_R2.fastq.gz"

  #Define the merged file
  FILE3="${MERGED_PATH}/${BASENAME}_filtered_aligned.merged.fastq"

  #Assembling with MEGAHIT
  megahit -1 "$FILE1" -2 "$FILE2" -r "$FILE3" -o "$OUTPUT_PATH/${BASENAME}" -t 12
done
