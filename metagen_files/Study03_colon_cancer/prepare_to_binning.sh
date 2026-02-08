#! /bin/bash

READS_DIR="./host_adapt_trimm_reads"
CONTIGS_DIR="./megahit_outputs"
OUT_PATH="bowtie_index_paired"

mkdir -p "$OUT_PATH"

for R1 in $READS_DIR/*.clean.pair.rmhost.1.fq.gz; do

    # Find mate R2
    R2="${R1/.clean.pair.rmhost.1.fq.gz/.clean.pair.rmhost.2.fq.gz}"

    # Name of sample (ex.: CRR224755)
    BASENAME=$(basename "$R1" .clean.pair.rmhost.1.fq.gz)

    # Directory where is the index of sample
    INDEX_DIR="$CONTIGS_DIR/$BASENAME"

    # Prefix of index
    INDEX_PREFIX="$INDEX_DIR/${BASENAME}_index"

    bowtie2 -x "$INDEX_PREFIX" -1 "$R1" -2 "$R2" -S "${BASENAME}.sam" -p 5
    samtools view -@ 4 -bS "${BASENAME}.sam" > "${BASENAME}.bam"
    samtools sort "${BASENAME}.bam" -o "${BASENAME}_sorted_index.bam"
    samtools index "${BASENAME}_sorted_index.bam"

    rm "${BASENAME}.bam" "${BASENAME}.sam"

    mv "${BASENAME}_sorted_index.bam" "$OUT_PATH"
    mv "${BASENAME}_sorted_index.bam.bai" "$OUT_PATH"

done
