#! /bin/bash

READS_DIR="./clean_reads_FR"
CONTIGS_DIR="./megahit_outputs"
OUT_PATH="bowtie_index_FR_paired"

mkdir -p "$OUT_PATH"

for R1 in $READS_DIR/*_filtered_aligned_R1.fastq.gz; do

    # Arquivo R2 correspondente
    R2="${R1/_filtered_aligned_R1.fastq.gz/_filtered_aligned_R2.fastq.gz}"

    # Nome da amostra (ex.: CRR224755)
    BASENAME=$(basename "$R1" _filtered_aligned_R1.fastq.gz)

    # Diretório onde está o índice dessa amostra
    INDEX_DIR="$CONTIGS_DIR/$BASENAME"

    # Prefixo do índice
    INDEX_PREFIX="$INDEX_DIR/${BASENAME}_index"

    bowtie2 -x "$INDEX_PREFIX" -1 "$R1" -2 "$R2" -S "${BASENAME}.sam"
    samtools view -bS "${BASENAME}.sam" > "${BASENAME}.bam"
    samtools sort "${BASENAME}.bam" -o "${BASENAME}_sorted_index.bam"
    samtools index "${BASENAME}_sorted_index.bam"

    rm "${BASENAME}.bam" "${BASENAME}.sam"

    mv "${BASENAME}_sorted_index.bam" "$OUT_PATH"
    mv "${BASENAME}_sorted_index.bam.bai" "$OUT_PATH"

done
