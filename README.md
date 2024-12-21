# Sialylation project of microbiome

![Bacteria%20sia](https://github.com/edtankian97/microbiota_sialylation/blob/teste/Bacteria%20sia.gif)

This project is related to my P.h.D.'s thesis which study is upon incorporation of sialic acid intestinal microbiome onto their cell wall. The main objective of the bioinformatic analysis is to
seek in proteomes of bacterias from NCBI that has a potential proteins linked to the process of sialic acid's incorporation. For this purpose, the bioinformatic part was divided in:

**1. Genome processing:** Retrieve information from NCBI, remotion of incomplete genomes and retrieve CheckM information.

**2. HMMER models:** Creation of proteins datasets and creation of models.

**3. Protein analysis:** Identification of sialylation pathway in filtered genomes.

**4. Downstream analysis:** Genomic analysis of bacterias' genomes that have sialylation with figures created.

## 1. Genome processing

### 1.1 Retrieve genome information

First of all, genomes were download with command wget at [NCBI](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt). 
 The original dataset is named **assembly_summary.txt** and be encountered in **genomes_download** directory.

### 1.2 Filtering NCBI retrieved dataset
```cd genomes_download
grep -c “Complete” assembly_summary.txt
grep "Complete" assembly_summary.txt > assembly_complete
cut -f1,8,9,20 assembly_complete > assembly_complete_summary.tsv #retrieve info that I want to get
```

### 1.3 Retrieve CheckM information
Please, check the ipynb file named as **Checkm_refseq_Reanalise_V2_R.ipynb**


