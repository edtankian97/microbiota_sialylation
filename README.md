# Sialylation project of microbiota

This project is related to my P.h.D.'s thesis which study is upon incorporation of sialic acid intestinal microbiota onto their cell wall. The main objective of the bioinformatic analysis is to
seek in proteomes of bacterias from NCBI that has a potential proteins linked to the process of sialic acid's incorporation. For this purpose, the bioinformatic part was divided in:

- 1. Genome processing: Retrieve information from NCBI, remotion of incomplete genomes and retrieve CheckM information.
- 2. HMMER models: Creation of proteins datasets and creation of models.
- 3. Protein analysis: Identification of sialylation pathway in filtered genomes.
- 4. Downstream analysis: Genomic analysis of bacterias' genomes that have sialylation with figures created.

## 1. Genome processing
First of all, genomes were download at NCBI site, which was retrieve with the following command: 
```console wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt ```.The original dataset is named as assembly_summary.txt and be encountered in genomes_download directory.
