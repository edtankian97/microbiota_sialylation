#! /bin/bash -v

#Define directories with HMM models and proteins sequences


#  extract their names
for mafft_file in ../Protein_database/CD_HIT/mafft_align/*; do
    # Extrair o nome do modelo (sem a extensão .fasta)
    mafft_name=$(basename "$mafft_file" .fasta)


        # Execute hmmbuild 
        output_file="${mafft_name}.hmm"
        hmmbuild  "$output_file" "$mafft_file" 
mv "$output_file" ../Protein_database/CD_HIT/mafft_align/hmm_models/
    done
