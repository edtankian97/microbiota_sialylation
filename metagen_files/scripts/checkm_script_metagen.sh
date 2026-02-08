#!/bin/bash

checkm lineage_wf \
  -x fa \
  -t 15 \
  -f lineage_log.txt \
  ./bins_paired/bins_paired_sialylation \
  checkm_result_bins_with_sia


checkm qa \
  ./checkm_result_bins_with_sia/lineage.ms \
  ./checkm_result_bins_with_sia \
  --tab_table \
  -f tabela_checkm.tsv