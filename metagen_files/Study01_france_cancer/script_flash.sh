#! /bin/bash/

# Dependencies: FLASH
# Install:
#   conda install -c bioconda flash
# Define paths

FASTQ_PATH="./clean_reads_FR"
OUTPUT_PATH="./outputs_flash"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_PATH"

# Loop through all *_R1.fastq.gz files in the FASTQ_PATH
for FILE1 in "$FASTQ_PATH"/*_R1.fastq.gz; do
  # Get the base name of the sample (remove the _R1.fastq.gz suffix)
  BASENAME=$(basename "$FILE1" _R1.fastq.gz)

  # Define the corresponding file for the _R2.fastq pair
  FILE2="${FASTQ_PATH}/${BASENAME}_R2.fastq.gz"

  # Check if the _R2.fastq file exists
  if [ -f "$FILE2" ]; then
    # Merging with FLASH
    flash -z -M 150 --to-stdout "$FILE1" "$FILE2" > "${OUTPUT_PATH}/${BASENAME}.merged.fastq"
    # -M 150: avoids the WARNING from FLASH saying that many sequences overlap more than the parameter, but it didn't change the merge percentage
  else
    echo "Pair for file $FILE1 not found: $FILE2" >&2
  fi
done
