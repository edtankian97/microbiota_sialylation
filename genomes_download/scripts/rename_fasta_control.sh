#!/bin/bash

# Itera sobre todos os arquivos no diretório com a extensão .faa
for arquivo in ../control_proteins/*.faa; do
    # Verifica se o arquivo existe e é um arquivo regular
    if [ -f "$arquivo" ]; then
        # Adiciona o nome do arquivo para cada linha que começa com '>'
        sed -i "s/^>/>$arquivo /" "$arquivo"
        echo "Nome do arquivo adicionado para cada linha no arquivo $arquivo"
    fi
done

