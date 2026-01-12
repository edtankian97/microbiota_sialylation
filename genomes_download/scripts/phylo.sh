#! bin/bash


phylophlan -i ../proteins/proteins_sialylation/final_complete_sialylation/ --db_type a -d phylophlan --databases_folder ../proteins/protein_tree/phylophlan_database/ --diversity medium --accurate -f ../proteins/protein_tree/genome_ed_config.cfg -o ../proteins/protein_tree/phylo_result/  --nproc 15  --verbose 2>&1 | tee ../proteins/protein_tree/logs/phylophlan_output_phylophlan.log
