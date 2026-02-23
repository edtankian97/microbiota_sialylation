#!/bin/bash

INPUT_FILE=./reads_sialylation/


for R1 in $INPUT_FILE/*_filtered_aligned_R1.fastq.gz
do
    SAMPLE=$(basename $R1 _filtered_aligned_R1.fastq.gz)

    R2=${R1/_R1/_R2}

    kraken2 \
        --db kraken2_DB \
        --paired $R1 $R2 \
        --threads 6 \
        --report ${SAMPLE}.report \
        --output ${SAMPLE}.kraken
done
