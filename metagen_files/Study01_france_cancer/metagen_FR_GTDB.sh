#! /bin/bash

export GTDBTK_DATA_PATH=../GTDB_data/split_package/release226/

gtdbtk classify_wf \
--extension fa \
  --genome_dir ./bins_paired/bins_paired_sialylation/ \
  --out_dir gtdbtk_out_FR_paired \
  --cpus 8
