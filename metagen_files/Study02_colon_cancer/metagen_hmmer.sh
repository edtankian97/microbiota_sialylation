#! /bin/bash

mkdir -p metagen_hmmer_results/

# Loop through the HMM models files with .hmm to extract their names
HMM_DIR="../../Protein_database/mafft_align/hmm_models/"

for hmm_file in "${HMM_DIR}"/*hmm; do
    # Extract models'name (remove .hmm)
    hmm_name=$(basename "$hmm_file" .hmm)

    # Loop through proteins files to extract their names
    for seq_file in ./bins_paired/all_prokka/*faa; do
        # Extract name of sequence (sem a extensão .faa)
        seq_name=$(basename "$seq_file" .faa)

        # Execute hmmsearch e redirect output to output file
        output_file="${hmm_name}_${seq_name}_output.tsv"
        hmmsearch --cpu 15 --domtblout "$output_file" "$hmm_file" "$seq_file"

 #Coverage calculation
        coverage_file="${output_file}_coverage"
        awk 'NR==4 {result = ((($19 - $18) + 1) / $3) * 100; printf "%s %s\n", result, ARGV[1]}' "$output_file" > "$coverage_file"

        # Add the content of line 12 as a new columm in coverage file
        paste -d' ' "$coverage_file" > temp_file && mv temp_file "$coverage_file"
        echo "Executado: $hmm_name contra $seq_name. Saída em $output_file"
        mv "$output_file"  ./metagen_hmmer_results/
 mv "$coverage_file" ./metagen_hmmer_results/

done
        done
