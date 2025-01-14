#! /bin/bash -v

#execute cd-hit and move resulted file

#for mixed CMP synthase
cd-hit -i ../Protein_database/CMP_synthase_mixed_database.fasta -o CMP_synthase_mixed_filter_100  -c 1.00 -n 5 
mv CMP_synthase_mixed_filter_100* ../Protein_database/CD_HIT/

#for review one
cd-hit -i ../Protein_database/CMP_synthase_review_database.fasta -o CMP_synthase_review_filter_100  -c 1.00 -n 5 
mv CMP_synthase_review_filter_100* ../Protein_database/CD_HIT/

# unreview CMP
cd-hit -i ../Protein_database/CMP_synthase_unreview_database.fasta -o CMP_synthase_unreview_filter_100  -c 1.00 -n 5 
mv CMP_synthase_unreview_filter_100* ../Protein_database/CD_HIT/

#sialiltransferase 

cd-hit -i ../Protein_database/sialiltransferase_mixed_database_old_gold.fasta -o old_gold_sialil_mixed_filter_100  -c 1.00 -n 5 
cd-hit -i ../Protein_database/sialiltransferase_review_database.fasta -o sialiltransferase_review_filter_100  -c 1.00 -n 5 
cd-hit -i ../Protein_database/sialiltransferase_unreview_database.fasta -o sialiltransferase_unreview_filter_100  -c 1.00 -n 5 
mv sialiltransferase* ../Protein_database/CD_HIT/


#polisialiltransferase 
cd-hit -i ../Protein_database/polisialil_database.fasta -o polisialil_filter_100  -c 1.00 -n 5 
mv polisialil_filter_100* ../Protein_database/CD_HIT/

#(Poly)O-acetiltransferase
cd-hit -i ../Protein_database/oacetil_plus_poli_mixed_database.fasta -o oacetil_plus_poli_mixed_filter_100  -c 1.00 -n 5 
mv oacetil_plus_poli_mixed_filter_100* ../Protein_database/CD_HIT/


#KpsM
cd-hit -i ../Protein_database/KpsM_mixed_database -o KpsM_mixed_filter_100  -c 1.00 -n 5 
mv KpsM_mixed_filter_100* ../Protein_database/CD_HIT/


#KpsT
cd-hit -i ../Protein_database/KpsM_mixed_database -o KpsM_mixed_filter_100  -c 1.00 -n 5 
mv KpsM_mixed_filter_100* ../Protein_database/CD_HIT/

#Rfah
cd-hit -i ../Protein_database/RfaH_mixed_database -o RfaH_mixed_filter_100  -c 1.00 -n 5 
mv RfaH_mixed_filter_100* ../Protein_database/CD_HIT/

