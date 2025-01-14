#! bin/bash

#Sequences alignment with mafft of sequence database

mafft --auto  --thread -5 ../Protein_database/CD_HIT/CMP_synthase_mixed_filter_100.fasta > CMP_synthase_mixed_mafft.fasta 
mv CMP_synthase_mixed_mafft.fasta  ../Protein_database/CD_HIT/mafft_align

mafft --auto  --thread -5 ../Protein_database/CD_HIT/CMP_synthase_unreview_filter_100.fasta > CMP_synthase_unreview_mafft.fasta 
mv CMP_synthase_unreview_mafft.fasta  ../Protein_database/CD_HIT/mafft_align

mafft --auto  --thread -5 ../Protein_database/CD_HIT/CMP_synthase_unreview_filter_100.fasta > CMP_synthase_unreview_mafft.fasta 
mv CMP_synthase_unreview_mafft.fasta  ../Protein_database/CD_HIT/mafft_align

mafft --auto  --thread -5 ../Protein_database/CD_HIT/old_gold_sialil_mixed_filter_100.fasta > old_gold_sialil_mixed_mafft.fasta 
mv old_gold_sialil_mixed_mafft.fasta   ../Protein_database/CD_HIT/mafft_align

mafft --auto  --thread -5 ../Protein_database/CD_HIT/sialiltransferase_review_filter_100.fasta > sialiltransferase_review_mafft.fasta 
mv sialiltransferase_review_mafft.fasta  ../Protein_database/CD_HIT/mafft_align

mafft --auto  --thread -5 ../Protein_database/CD_HIT/sialiltransferase_unreview_filter_100.fasta > sialiltransferase_unreview_mafft.fasta 
mv sialiltransferase_unreview_mafft.fasta  ../Protein_database/CD_HIT/mafft_align

mafft --auto  --thread -5 ../Protein_database/CD_HIT/polisialil_filter_100.fasta > polisialil_mafft.fasta 
mv polisialil_mafft.fasta  ../Protein_database/CD_HIT/mafft_align

mafft --auto  --thread -5 ../Protein_database/CD_HIT/oacetil_plus_poli_mixed_filter_100.fasta > oacetil_plus_poli_mixed_mafft.fasta 
mv oacetil_plus_poli_mixed_mafft.fasta  ../Protein_database/CD_HIT/mafft_align

mafft --auto  --thread -5 ../Protein_database/CD_HIT/KpsM_mixed_filter_100.fasta > KpsM_mixed_mafft.fasta 
mv KpsM_mixed_mafft.fasta  ../Protein_database/CD_HIT/mafft_align

mafft --auto  --thread -5 ../Protein_database/CD_HIT/KpsT_mixed_filter_100.fasta > KpsT_mixed_mafft.fasta 
mv KpsT_mixed_mafft.fasta  ../Protein_database/CD_HIT/mafft_align

mafft --auto  --thread -5 ../Protein_database/CD_HIT/RfaH_mixed_filter_100.fasta > RfaH_mixed_mafft.fasta 
mv RfaH_mixed_mafft.fasta  ../Protein_database/CD_HIT/mafft_align

