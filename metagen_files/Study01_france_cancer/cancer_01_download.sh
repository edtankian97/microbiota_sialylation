#! /bin/bash -v

while read accession; do
        fasterq-dump --split-files "${accession}" --outdir france_fastq_reads/
    done < all_cancer_download.sh
