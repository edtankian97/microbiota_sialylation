#! /bin/bash


# Diretório onde estão os diretórios das amostras
SAMPLES_DIR="./megahit_outputs/"

# Percorre cada subdiretório
for DIR in "$SAMPLES_DIR"/*; do
    # Verifica se é diretório
    if [ -d "$DIR" ]; then

        # Encontra o arquivo de contigs
        CONTIGS=$(find "$DIR" -maxdepth 1 -type f -name "*.fa")


        BASENAME=$(basename "$DIR")
        INDEX_PREFIX="$DIR/${BASENAME}_index"

        bowtie2-build "$CONTIGS" "$INDEX_PREFIX"
    fi
done

echo "Finalizado!"
