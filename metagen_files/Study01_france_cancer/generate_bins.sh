#! /bin/bash

# Directories
BAM_DIR="bowtie_index_FR_paired"
CONTIG_DIR="./megahit_outputs/"
OUT_DIR="bins_paired"

mkdir -p "$OUT_DIR"

# Loop for each BAM
for BAM in "$BAM_DIR"/*_sorted_index.bam; do

    # Extract name of sample
    SAMPLE=$(basename "$BAM" _sorted_index.bam)


 CONTIG="$CONTIG_DIR/$SAMPLE/final.contigs.fa"

    # Generate file of metabat binning depth
    DEPTH="$OUT_DIR/depth_${SAMPLE}.txt"
    jgi_summarize_bam_contig_depths \
        --outputDepth "$DEPTH" \
        "$BAM"

    # Create directory for bins
    BIN_OUT="$OUT_DIR/$SAMPLE"
    mkdir -p "$BIN_OUT"

    # Run metabat2
    metabat2 \
        -i "$CONTIG" \
        -a "$DEPTH" \
        -o "$BIN_OUT/bin"

done
