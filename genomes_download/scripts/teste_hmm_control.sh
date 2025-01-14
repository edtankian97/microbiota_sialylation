#! /bin/bash -v

#Define directories with HMM models and proteins sequences 

seq_dir=../control_proteins/

# Loop through the HMM models files with .hmm to extract their names
for hmm_file in ../Protein_database/CD_HIT/mafft_align/hmm_models/*; do
    # Extrair o nome do modelo (sem a extensão .hmm)
    hmm_name=$(basename "$hmm_file" .hmm)

    # Loop through proteins files to extract their names
    for seq_file in "$seq_dir"/*.faa; do
        # Extrair o nome do arquivo de sequência (sem a extensão .faa)
        seq_name=$(basename "$seq_file" .faa)

        # Execute hmmsearch e redirect a saída para um arquivo de saída
        output_file="${hmm_name}_${seq_name}_output.tsv"
        hmmsearch --cpu 15 --domtblout "$output_file" "$hmm_file" "$seq_file" 

        #Coverage calculation
        coverage_file="${output_file}_coverage"
        awk 'NR==4 {result = ((($19 - $18) + 1) / $3) * 100; printf "%s %s\n", result, ARGV[1]}' "$output_file" > "$coverage_file"
	
        # Add the content of line 12 as a new columm in coverage file
        paste -d' ' "$coverage_file" > temp_file && mv temp_file "$coverage_file" 
        echo "Executado: $hmm_name contra $seq_name. Saída em $output_file"
	mv "$output_file" ../control_proteins/HMMER_CONTROL_RESULTS
	mv "$coverage_file" ../control_proteins/HMMER_CONTROL_RESULTS
    done
done


