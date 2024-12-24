# Sialylation project of microbiome

![Bacteria%20sia](https://github.com/edtankian97/microbiota_sialylation/blob/teste/Bacteria%20sia.gif)

This project is related to my P.h.D.'s thesis which study is upon incorporation of sialic acid intestinal microbiome onto their cell wall. The main objective of the bioinformatic analysis is to
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
pwd #you should be in the directory **genomes_download**
sed '1d' remain_CheckM_data_complete_with_NA.tsv > remain_CheckM_data_complete_with_NA_ID.tsv
datasets download genome accession --inputfile remain_CheckM_data_complete_with_NA_ID.tsv  #Download genomes filesos .fna
unzip ncbi_dataset.zip && mv ncbi_dataset/ remain_CheckM/
rm ncbi_dataset.zip
find ./remain_CheckM/data/GCF_000*/ -type f -iname "*.fna" -exec mv -v "{}" ./remain_CheckM/ \;
checkm2 predict --threads 5 --input  remain_CheckM/ --output-directory checkm2_result
```
Now return to the script **Checkm_refseq_Reanalise_V2_R.ipynb**

## 2 

## 3. Protein analysis

### 3.1 Protein analysis: Control proteomes

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
