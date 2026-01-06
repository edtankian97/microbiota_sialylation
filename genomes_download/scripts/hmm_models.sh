#! /bin/bash -v

#Define directories with HMM models and proteins sequences

OUT_PATH=../Protein_database/mafft_align/hmm_models/

mkdir -p $OUT_PATH

#  extract their names
for mafft_file in ../Protein_database/mafft_align/*; do
    # Extract name of file (without extension .fasta)
    mafft_name=$(basename "$mafft_file" _mafft.fasta)


        # Execute hmmbuild 
        output_file="${mafft_name}.hmm"
        hmmbuild  "$output_file" "$mafft_file" 
mv "$output_file" $OUT_PATH
    done
