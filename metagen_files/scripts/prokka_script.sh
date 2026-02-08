#! /bin/bash

export PROKKA_NO_BLAST_CHECK=1

for f in ./bins_paired/ERR*/*.fa; do
    # extrai só o nome do arquivo sem caminho e sem extensão
    base=$(basename "$f" .fa)
    dir=$(dirname "$f")
    prokka "$f" \
      --outdir "${dir}/${base}.prokka" \
     --prefix "PROKKA_${base}" \
     --force

    echo "Finalizado: $f"
done
