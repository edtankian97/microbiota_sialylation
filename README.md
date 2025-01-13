# Sialylation project of microbiota

![Bacteria%20sia](https://github.com/edtankian97/microbiota_sialylation/blob/teste/Bacteria%20sia.gif)

This project is related to my P.h.D.'s thesis which study is upon incorporation of sialic acid intestinal microbiota onto their cell wall. The main objective of the bioinformatic analysis is to
seek in proteomes of bacterias from NCBI that has a potential proteins linked to the process of sialic acid's incorporation. For this purpose, the bioinformatic part was divided in:

**1. Genome processing:** Retrieve information from NCBI, remotion of incomplete genomes and retrieve CheckM information.

**2. HMMER models:** Creation of proteins datasets and creation of models.

**3. Protein analysis:** Identification of sialylation pathway in filtered genomes/proteomes and control proteomes.

**4. Downstream analysis:** Genomic analysis of bacterias' genomes that have sialylation with figures created.

## 1. Genome processing

### 1.1 Retrieve genome information

First of all, genomes were download with command wget at [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt). 
 The original dataset is named **assembly_summary.txt** and be encountered in **genomes_download** directory.

### 1.2 Filtering NCBI retrieved dataset
```
cd genomes_download
grep -c “Complete” assembly_summary.txt
grep "Complete" assembly_summary.txt > assembly_complete
cut -f1,8,9,20 assembly_complete > assembly_complete_summary.tsv #retrieve info that I want to get
```

### 1.3 Retrieve CheckM information
First of all, start with the ipynb file named as **Checkm_refseq_Reanalise_V2_R.ipynb** which will be in this path=genomes_download/scripts/jupyter_scripts/
If you do not have installed miniconda or anaconda yet, please follow the instructions in this [link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).

**-Creation of conda environment and installation of ncbi_datasets. More information about ncbi_datasets, click this [link](https://github.com/ncbi/datasets)**
```
conda create -n ncbi_datasets python=3.8 #creation of the anaconda environment: Digit y or yes to continue the installation. If error occurs, you might update the python version.
conda activate ncbi_datasets #Activation of the environment. Do this after creation of the environment
conda install -c conda-forge ncbi-datasets-cli
```
**-retrieve missing data of completeness from ncbi_datasets**
```
xargs -a GCF_complete_without_checkM.txt -I {} datasets summary genome accession {} --as-json-lines | dataformat tsv genome --fields organism-name,accession,checkm-completeness,checkm-contamination > remain_CheckM_data.tsv
mv remain_CheckM_data.tsv remain_CheckM_data_complete.tsv
```
after this, continue the Part 02 of script: **Checkm_refseq_Reanalise_V2_R.ipynb**

**Run CheckM2 to get completeness and contamination of missing data. CheckM2 uses artificial language to predict completeness of genomes. 
```
conda activate ncbi_datasets
conda install -c bioconda -c conda-forge checkm2 
checkm2 -h
checkm2 database --download #download checkm2's standart database
```
**Download missing genomes for contamination analysis and run CheckM2**
```
pwd #you should be in the directory genomes_download
sed '1d' remain_CheckM_data_complete_with_NA.tsv > remain_CheckM_data_complete_with_NA_ID.tsv
datasets download genome accession --inputfile remain_CheckM_data_complete_with_NA_ID.tsv  #Download genomes files ".fna"
unzip ncbi_dataset.zip && mv ncbi_dataset/ remain_CheckM/
rm ncbi_dataset.zip
find ./remain_CheckM/data/GCF_000*/ -type f -iname "*.fna" -exec mv -v "{}" ./remain_CheckM/ \;
checkm2 predict --threads 5 --input  remain_CheckM/ --output-directory checkm2_result
```
Now return to the script **Checkm_refseq_Reanalise_V2_R.ipynb**

## 2 HMMER models

With all protein fasta downloaded, remove duplicated sequences with CD-HIT, which can be downloaded following this [guide](https://github.com/weizhongli/cdhit)
After installation, check if the PATH for CD-HIT executable is recognize in your system. If not, you might add PATH to CD-HIT executable folder.
(Example usage of cd-hit: cd-hit -i [PROTEIN_FASTA_NAME] -o [CD_HIT_ENZYME_NAME_MODE_TYPE_OUTPUT_FILE]  -c 1.00 -n 5 ).
```
cd Protein_database
mkdir CD_HIT
bash ../scripts/CD_HIT_script.sh
```
After this, it's turn to do an alignment. For this purpose, follow [mafft](https://mafft.cbrc.jp/alignment/software/) installation.  
(Example of mafft usage: mafft --auto [CD_HIT_ENZYME_NAME_MODE_TYPE_OUTPUT_FILE] > [ENZYME_NAME_MODE_TYPE_OUTPUT_FILE]_mafft.fasta).
```
cd CD_HIT
mkdir mafft_align
bash ../../scripts/mafft_align.sh 
```
In the end, let's construct protein models with [HMMER](https://github.com/EddyRivasLab/hmmer).In this link aside, I've installed the uppermost recent version #3.4 by downloading the raw file.
```
cd mafft_align
mkdir hmm_models
bash ../../../scrips/hmm_models.sh #results will be located in hmm_models folder
cd ../../../../ #you must be located in genomes_download
```
## 3. Protein analysis

### 3.1 Protein analysis: Control proteomes

First, we will download some proteomes based on NCBI ID which are presented in **control_proteomes.txt**. Other ones are already available in **control_proteins** folder
```
conda activate ncbi_datasets
cd ./control_proteins/ 
datasets download genome accession --inputfile control_proteomes.txt --include protein --filename control.zip
unzip control.zip -d control
```
After download, process the files
```
bash ../../scripts/rename_control_files.sh #rename the files based on their directories
find ./control/ncbi_dataset/data/GCF*/ -type f -iname "*.faa" -exec mv -v "{}" ./ \; #move files
while read line; do eval mv $line; done < files.txt #rename with species names
bash ../../scripts/rename_fasta_control.sh #rename fasta header
```
Create directories to organize the results.
```
mkdir HMMER_CONTROL_RESULTS && cd HMMER_CONTROL_RESULTS
mkdir CMP sialil polisia kpsM_T Rfah
mkdir CMP/mixed CMP/review CMP/unreview sialil/mixed sialil/review sialil/unreview 
bash ../../scripts/teste_hmm_control.sh
```
Join all output files for each enzyme
```
cat CMP_mixed*_output.tsv > all_CMP_mixed_control_output.tsv
cat CMP_review*_output.tsv > all_CMP_review_control_output.tsv
cat CMP_unreview*_output.tsv > all_CMP_unreview_control_output.tsv
cat sialil_review*_output.tsv > all_sialil_review_control_output.tsv
cat sialil_unreview*_output.tsv > all_sialil_unreview_control_output.tsv
cat all_sialil_old_one*_output.tsv > all_sialil_mixed_control_output.tsv
cat polisialil*_output.tsv > all_polisialil_control_output.tsv
cat Rfah_1*_output.tsv > all_Rfah_1_control_output.tsv
cat kpsM_mafft*_output.tsv > all_kpsM_control_output.tsv
cat kpsT_mafft*_output.tsv > all_kpsT_control_output.tsv
```
Do the same thing for coverage files
```
cat CMP_mixed*_output.tsv_coverage > all_CMP_mixed_control_output_coverage.tsv
cat CMP_review*_output.tsv_coverage > all_CMP_review_control_output_coverage.tsv
cat CMP_unreview*_output.tsv_coverage > all_CMP_unreview_control_output_coverage.tsv
cat sialil_review*_output.tsv_coverage > all_sialil_review_control_output_coverage.tsv
cat sialil_unreview*_output.tsv_coverage > all_sialil_unreview_control_output_coverage.tsv
cat all_sialil_old_one*_output.tsv_coverage > all_sialil_mixed_control_output_coverage.tsv
cat polisialil*_output.tsv_coverage > all_polisialil_control_output_coverage.tsv
cat Rfah_1*_output.tsv_coverage > all_Rfah_1_control_output_coverage.tsv
cat kpsM_mafft*_output.tsv_coverage > all_kpsM_control_output_coverage.tsv
cat kpsT_mafft*_output.tsv_coverage > all_kpsT_control_output_coverage.tsv
```
Process output file
```
sed -i '/#/d' *_output.tsv
sed -i 's/ \{1,\}/\t/g' *_output.tsv 
```
Move the results to their respectives directories
```
mv CMP_mixed* CMP/mixed
mv CMP_review* CMP/review
mv CMP_unreview* CMP/unreview
mv sialiltransferase_review* sialil/review
mv sialiltransferase_unreview* sialil/unreview
mv all_sialil_old_one* sialil/mixed
mv polisialil* polisia/
mv Rfah_1* Rfah/
mv kpsM_mafft* kpsM_T/
mv kpsT_mafft* kpsM_T/
```
### 3.2 Protein analysis: NCBI analysis
After filtration with CheckM, the proteomes were downloaded from NCBI
```
cut -f2 checkm_filter_v2_complete.tsv > checkm_filter_v2_complete_ID.tsv
sed -i '1d' checkm_filter_v2_complete_ID.tsv
datasets download genome accession --inputfile checkm_filter_v2_complete_ID.tsv --include protein --dehydrated  --filename proteins.zip
unzip proteins.zip -d proteins
datasets rehydrate --directory proteins
```
After the download, proteins files were renamed with their own directory name 
```
cd scripts/
bash rename_files.sh
```
With the proteins with their respective names, you can move them to the **proteins** directory
```
cd .. # you should be in genomes_download directory
find proteins/ncbi_dataset/data/GCF*/ -type f -iname "*.faa" -exec mv -v "{}" ./proteins/ \;
```
Proteomes are now in **proteins** directory and then you can edit their fasta header, so we can identify them later on HMMER analysis.
For this, do the following command:
```
cd scripts/
bash rename_fasta.sh
```
Now the proteins files are ready to be analised. Now, let's organize subdiretories to store the results for each enzyme of sialylation process
```
cd ..
mkdir HMMER_analysis
cd HMMER_analysis
mkdir CMP_synthase sialiltrans polisialiltrans o_acetiltrans_plus_poli kpsM kpsT RfaH
cd ..
```
Scripts for each enzyme model will be available with their respective names in **scripts** directory

**CMP_synthase**
```
cd scripts/
bash CMP_hmm.sh
```
**Sialiltransferase**
```
bash Sialiltrans_hmm.sh
```
**Polisialiltransferase**
```
bash polisialiltrans_hmm.sh
```
**(Poly)O-acetiltransferase**
```
bash o_acetiltrans_poli_hmm.sh
```
**KpsM**
```
bash KpsM_hmm.sh
```
**KpsT**
```
bash KpsT_hmm.sh
```
**RfaH**
```
bash RfaH_hmm.sh
```

With all done, now it's time to start the process of coverage, e-value and bit-score filter
```
cd ..
cd HMMER_analysis/
find ./CMP_synthase/ -type f -name '*coverage*' -exec cat {} + > CMP_coverage.tsv
find ./sialiltrans/ -type f -name '*coverage*' -exec cat {} + > sialil_coverage.tsv
find ./polisialiltrans/ -type f -name '*coverage*' -exec cat {} + > polisialil_coverage.tsv
find ./o_acetiltrans_plus_poli/ -type f -name '*coverage*' -exec cat {} + > o_acetiltrans_plus_poli_coverage.tsv
find ./kpsM/ -type f -name '*coverage*' -exec cat {} + > kpsM_coverage.tsv
find ./kpsT/ -type f -name '*coverage*' -exec cat {} + > kpsT_coverage.tsv
find ./RfaH/ -type f -name '*coverage*' -exec cat {} + > RfaH_coverage.tsv
```
Now, go to the script **hmm_process.ipynb** which is loccated in the path: microbial_sialylation/genomes_download/scripts/jupyter_scripts/ and follow the script.

After part 01 with coverage assessment, follow this for each core enzyme: CMP synthase, sialiltransferase and polisialiltransferase. 
**CMP synthase**
```
cd ./HMMER_analysis/CMP_synthase/
sed -i '1d' CMP_complete_cover_filter_40.tsv
cut -f2 CMP_complete_cover_filter_40.tsv
mkdir filter_cover_CMP
for file in $(cat ./CMP_complete_cover_filter_40.tsv); do mv "$file" ./filter_cover_CMP/; done
cd filter_cover_CMP
#output tsv
find . -type f -name '*protein_output*' -exec cat {} + > new_file.tsv
#process output file
sed -i '/#/d' new_file.tsv
sed  's/ \{1,\}/\t/g' new_file.tsv > file_output_CMP.tsv
cd ../../../
```
**Sialiltransferase**
```
cd ./HMMER_analysis/sialil/
sed -i '1d' sialil_complete_cover_filter_40.tsv
cut -f2 sialil_complete_cover_filter_40.tsv
mkdir filter_cover_sialil
for file in $(cat ./sialil_complete_cover_filter_40.tsv); do mv "$file" ./filter_cover_sialil/; done
cd filter_cover_sialil
#output tsv
find . -type f -name '*protein_output*' -exec cat {} + > new_file.tsv
#process output file
sed -i '/#/d' new_file.tsv
sed  's/ \{1,\}/\t/g' new_file.tsv > file_output_sialil.tsv
```
**polisialiltransferase**
```
cd ./HMMER_analysis/polisiatrans/
sed -i '1d' polisia_complete_cover_filter_40.tsv
cut -f2 polisia_complete_cover_filter_40.tsv
mkdir filter_cover_polisia
for file in $(cat ./polisia_complete_cover_filter_40.tsv); do mv "$file" ./filter_cover_polisia/; done
cd filter_cover_polisia
#output tsv
find . -type f -name '*protein_output*' -exec cat {} + > new_file.tsv
#process output file
sed -i '/#/d' new_file.tsv
sed  's/ \{1,\}/\t/g' new_file.tsv > file_output_polisia.tsv
```

# 4. Downstream analysis

## 4.1 Datasets for plots
This topic and subtopics forwards are about how to get data that will be important to create the plots.
### 4.1.1 Get information of complete genomes
```
#Retrieve info of complete genomes of geo_loc, isolation_source and other info
cd ../../../ #you must be located in genomes_download
cut -f3 checkm_filter_v2_complete.tsv > checkm_filter_v2_complete_ID.tsv
sed -i '1d' checkm_filter_v2_complete_ID.tsv
datasets download genome accession --inputfile checkm_filter_v2_complete_ID.tsv --dehydrated
unzip -v ncbi_dataset.zip
dataformat tsv genome --package ncbi_dataset.zip --fields accession,assminfo-biosample-geo-loc-name,assminfo-biosample-host,assminfo-biosample-host-disease,assminfo-biosample-isolation-source,assminfo-biosample-source-type > accession_list.tsv
dataformat tsv genome --package ncbi_dataset.zip --fields accession,assminfo-biosample-geo-loc-name,assminfo-biosample-host,assminfo-biosample-host-disease,assminfo-biosample-isolation-source,assmstats-gc-percent,assmstats-total-sequence-len,organelle-assembly-name,organism-name,organism-tax-id,annotinfo-featcount-gene-protein-coding > accession_complete_fields.tsv
```

### 4.1.2 Information of genomes with sialylation pathway

**After hmm_process analysis of core enzymes, do the following to get information of genomes that have sialylation pathway**
```
#join files
cat comm_CMP_sia_poli.tsv not_comm_CMP_sia_poli.tsv not_commm_CMP_sia_poli_1.tsv > all_commom.tsv
#After join the files, remove header "X1"
fgrep -v X1 all_commom.tsv > all_commom_1.tsv
#remove _protein.faa
sed 's/_protein.faa//' all_commom_1.tsv > all_commom_1_modified.tsv
#get info
datasets download genome accession --inputfile all_commom_1_modified.tsv --dehydrated
dataformat tsv genome --package ncbi_dataset.zip  > comm_complete_genomes_dataset_tbl.tsv
#select desired fields
dataformat tsv genome --package ncbi_dataset.zip --fields accession,assminfo-biosample-geo-loc-name,assminfo-biosample-host,assminfo-biosample-host-disease,assminfo-biosample-source-type,assmstats-gc-percent,assmstats-total-sequence-len,organelle-assembly-name,organism-name,organism-tax-id > comm_complete_genomes_dataset_fields.tsv
#move to the folder

```
### 4.1.3 Taxonomy information
```
#take desired columm
cut -f10 comm_complete_genomes_dataset_fields.tsv > comm_sia_genomes_tax_id
sed -i '1d' comm_sia_genomes_tax_id #remove header
#retrieve taxonomy information
datasets download taxonomy taxon --inputfile comm_sia_genomes_tax_id --filename taxonomy.zip
unzip taxonomy.zip
mv taxonomy/ncbi_dataset/data/taxonomy_summary.tsv ./plots_data/
```
### 4.1.4 Phylogenetic tree

Fist of all, we will set up unique species' proteomes for the phylogenetic tree. First open R or Rstudio and copy and paste the following code.

```r
library(readr)
library(dplyr)
library(stringr)


getwd()
select_1 <- read_tsv("L:/comm_complete_genomes_dataset_fields.tsv")
head(select_1)

select_2 <- select_1 %>% select(`Assembly Accession`, `Organism Name`)
select_2$extract_name <- word(select_2$`Organism Name`,1,2)

#select only rows with new species ("sp.")
select_sp <- select_2[grep("sp.$", select_2$extract_name), ]

#select only rows without new species ("sp.")
no_sp <- select_2[-grep("sp.$", select_2$extract_name), ]

#see uniques species in extract_name collumn
select_3 <- distinct(no_sp,extract_name, .keep_all = TRUE)

all_select <- rbind(select_3,select_sp)

select_ID <- all_select%>% select(`Assembly Accession`)

write_tsv(select_ID, "microbiota_sialylation/genomes_download/proteins/proteins_unique_ID.tsv")


Select representative species 
```
cd ./proteins/
sed -i '1d' proteins_unique_ID.tsv
mkdir proteins_comm_sia/ && mkdir proteins_comm_sia/proteins_unique_comm_sia
for file in $(cat ./proteins_unique_ID.tsv); do cp "$file" ./proteins_comm_sia/proteins_unique_comm_sia/; done
cd ./proteins_comm_sia/
ls ./proteins_unique_comm_sia > representative_species.txt
sed -i '1,2d' representative_species.txt
head -n -1  representative_species.txt > representative_species_v1.txt
sed 's/_protein.faa//g' representative_species_v1.txt > representative_species_modified_v6.txt
sed 's/$/*/g' representative_species_modified_v6.txt > representative_species_modified_v5.txt
sed -i 's/^/*/g' representative_species_modified_v5.txt
cp representative_species_modified_v6.txt ../../genomes_download/plots_data/itol/
```

To generate a tree from phylophlan, first you must download it. Check this [link](https://github.com/biobakery/phylophlan) with the procedures.
```
conda create -n "phylophlan" -c bioconda phylophlan=3.1.1
conda activate phylophlan
```
**phylophlan database setup and installation**
I followed instructions upon this [link](https://github.com/biobakery/phylophlan/wiki#databases). I followed the option 2 and installed the **phylophlan** database. Download **phylophlan_databases.txt** and follow the instructions below.
```
cd proteins/
mkdir protein_tree && cd protein_tree
mkdir phylophlan_database && cd phylophlan_database
cat phylophlan_databases.txt # copy and paste one of the links to **phylophlan** database (not amphora)
tar -xf phylophlan.tar
bunzip2 -k phylophlan/phylophlan.bz2
cd ..
```
**phylophlan configuration file**
```
phylophlan_write_config_file.py --db_aa diamond --map_aa diamond --msa mafft \
--trim trimal --tree1 fasttree --tree2 raxml -o genome_ed_config.cfg \
--db_type a
```
**Run script of phylophlan analysis**
```
cd ../../scripts/
bash phylo.sh #make sure phylophlan's conda environment is activated
```
phylophlan generates a lot of files, but the most important is refine tree called **RAxML_result.proteins_unique_comm_sia_refined.tre** which was used for tree annotation with iTOL.
This output is present in **genomes_download** folder

## 4.2 Genome information

Follow the script **retrieve_genome_info.ipynb** which is loccated in the path: microbial_sialylation/genomes_download/scripts/jupyter_scripts/

## 4.3 Host distribution

Follow the script **host_distribution.ipynb** which is loccated in the path: microbial_sialylation/genomes_download/scripts/jupyter_scripts/

## 4.4 Species distribution

Follow the script **pie_data.ipynb** which is loccated in the path: microbial_sialylation/genomes_download/scripts/jupyter_scripts/

## 4.5 iTOL annotation

Follow the **Part 02** of hmm_process.ipynb script.
For the other proteins (KpsT, KpsM and RfaH), I followed in another way: subsets for representative species of the dataset. This will be useful for iTOL annotation: 
**4.5 iTOL annotation**

**KpsM**
```
cp representative_species_modified_v5.txt  HMMER_analysis/kpsM/
cd HMMER_analysis/kpsM/
mkdir representative_species_kpsM
for file in $(cat ./representative_species_modified_v5.txt); do mv "$file" ./representative_species_kpsM/; done
cd representative_species_kpsM/
find ./ -type f -name '*coverage*' -exec cat {} + > kpsM_coverage.tsv
```
Process resulted coverage file in hmm_process.ipynb script and then continue
```
#remove first line and first columm to get only ID
sed -i '1d' kpsM_representative_cover_filter_40.tsv
cut -f2 kpsM_representative_cover_filter_40.tsv > kpsM_representative_cover_filter_40_ID.tsv 
mkdir filtered_kpsM
#transfer files to another folder
for file in $(cat ./kpsM_representative_cover_filter_40_ID.tsv); do mv "$file" ./filtered_kpsM/; done
cd filtered_kpsM
find . -type f -name '*protein_output*' -exec cat {} + > all_kpsM_protein.tsv
#process output file
sed -i '/#/d' all_kpsM_protein.tsv
sed -i 's/ \{1,\}/\t/g' all_kpsM_protein.tsv 
```
Return to the script again. Final resulted file will be located in microbiota_sialylation/genomes_download/plots_data/itol/

**KpsT**
```
cp representative_species_modified_v5.txt  HMMER_analysis/kpsT/
cd HMMER_analysis/kpsT/
mkdir representative_species_kpsT
for file in $(cat ./representative_species_modified_v5.txt); do mv "$file" ./representative_species_kpsT/; done
cd representative_species_kpsT/
find ./ -type f -name '*coverage*' -exec cat {} + > kpsT_coverage.tsv
```
Process resulted coverage file in hmm_process.ipynb script and then continue
```
#remove first line and first columm to get only ID
sed -i '1d' kpsT_representative_cover_filter_40.tsv
cut -f2 kpsT_representative_cover_filter_40.tsv > kpsT_representative_cover_filter_40_ID.tsv 
mkdir filtered_kpsT
#transfer files to another folder
for file in $(cat ./kpsT_representative_cover_filter_40_ID.tsv); do mv "$file" ./filtered_kpsT/; done
cd filtered_kpsT
find . -type f -name '*protein_output*' -exec cat {} + > all_kpsT_protein.tsv
#process output file
sed -i '/#/d' all_kpsT_protein.tsv
sed -i 's/ \{1,\}/\t/g' all_kpsT_protein.tsv 
```
Return to the script again. Final resulted file will be located in microbiota_sialylation/genomes_download/plots_data/itol/

**(Poly)O-acetiltransferase**
```
cp representative_species_modified_v5.txt  HMMER_analysis/o_acetiltrans_plus_poli/
cd HMMER_analysis/o_acetiltrans_plus_poli/
mkdir representative_oacetil
for file in $(cat ./representative_species_modified_v5.txt); do mv "$file" ./representative_oacetil/; done
cd representative_oacetil/
find ./ -type f -name '*coverage*' -exec cat {} + > all_oacetil_cover
```
Process resulted coverage file in hmm_process.ipynb script and then continue
```
#remove first line and first columm to get only ID
sed -i '1d' oacetil_representative_cover_filter_40.tsv
cut -f2 oacetil_representative_cover_filter_40.tsv > oacetil_representative_cover_filter_40_ID.tsv 
mkdir filtered_oacetil
#transfer files to another folder
for file in $(cat ./oacetil_representative_cover_filter_40_ID.tsv); do mv "$file" ./filtered_oacetil/; done
cd filtered_oacetil
find . -type f -name '*protein_output*' -exec cat {} + > all_oacetil_plus_poli.tsv
#process output file
sed -i '/#/d' all_oacetil_plus_poli.tsv
sed -i 's/ \{1,\}/\t/g' all_oacetil_plus_poli.tsv 
```
Return to the script again. Final resulted file will be located in microbiota_sialylation/genomes_download/plots_data/itol/

**RfaH**
```
cp representative_species_modified_v5.txt  HMMER_analysis/RfaH/
cd HMMER_analysis/RfaH/
mkdir representative_RfaH
for file in $(cat ./representative_species_modified_v5.txt); do mv "$file" ./representative_RfaH/; done
cd representative_RfaH/
find ./ -type f -name '*coverage*' -exec cat {} + > rfah_coverage
```
Process resulted coverage file in hmm_process.ipynb script and then continue
```
#remove first line and first columm to get only ID
sed -i '1d' rfah_representative_cover_filter_40.tsv
cut -f2 rfah_representative_cover_filter_40.tsv > rfah_representative_cover_filter_40_ID.tsv 
mkdir filter_rfaH
#transfer files to another folder
for file in $(cat ./rfah_representative_cover_filter_40_ID.tsv); do mv "$file" ./filter_rfaH/; done
cd filter_rfaH
find . -type f -name '*protein_output*' -exec cat {} + > all_rfah_tsv_w_o_cruz.tsv
#process output file
sed -i '/#/d' all_rfah_tsv_w_o_cruz.tsv
sed -i 's/ \{1,\}/\t/g' all_rfah_tsv_w_o_cruz.tsv 
```
Return to the script again. Final resulted file will be located in microbiota_sialylation/genomes_download/plots_data/itol/


Follow the script **itol_notation.ipynb** which is loccated in the path: microbial_sialylation/genomes_download/scripts/jupyter_scripts/
After this, each dataset created in the script was concatenated with iTol dataset
```
#remove header
sed -i '1d' kpsT_represent_itol_EANS.tsv
sed -i '1d' oacetil_itol_represent_EANS.tsv
sed -i '1d' kpsM_itol_represent_EANS.tsv
sed -i '1d' vfdb_itol_represent_EANS.tsv
sed -i '1d' representative_labels_itol.tsv
```
#Ensure to move the to plots_data/itol folder, if is not there. 
There will be files that are from iTOL called 'datasets' in general
```
 #Join files
cat kpst_itol_representative.txt kpsT_represent_itol_EANS.tsv > final_kpsT_itol.txt
cat oacetil_itol_representative.txt oacetil_itol_represent_EANS.tsv > final_oacetil_itol.txt
cat kpsM_itol_representative.txt kpsM_itol_represent_EANS.tsv > final_kpsM_itol.txt
cat vfdb_itol_representative_EANS.txt vfdb_itol_represent_EANS.tsv > final_vfdb_itol.txt
cat ed_tree_label.txt representative_labels_itol.tsv > ed_tree_label_represent.txt
```
Now you can upload the files in iToL site.
