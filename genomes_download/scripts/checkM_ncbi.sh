#!/bin/bash

checkm lineage_wf \
  -t 6 \
  -f lineage_log.txt \
  ./remain_CheckM/ \
  ../checkm_result_ncbi/

checkm qa \
  ../checkm_result_ncbi/lineage.ms \
  ../checkm_result_ncbi/ \
  --tab_table \
  -f tabela_checkm_ncbi.tsv