#! /bin/bash



for f in ./ERR*/*.fa; do
    # extrai só o nome do arquivo sem caminho e sem extensão
    base=$(basename "$f" .fa)

    prokka "$f" \
        --outdir "${f}/${base}.prokka" \
        --prefix "PROKKA_${base}"

    echo "Finalizado: $f"
done
