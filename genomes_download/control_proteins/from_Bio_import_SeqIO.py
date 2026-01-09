from Bio import SeqIO

# Path to your GenBank file
gbk_file = "Escherichia_coli_ATCC_23511.gbk"

# List to store all protein sequences
proteins = []

# Parse GenBank file
for record in SeqIO.parse(gbk_file, "genbank"):
    for feature in record.features:
        if feature.type == "CDS" and "translation" in feature.qualifiers:
            protein_seq = feature.qualifiers["translation"][0]
            proteins.append(protein_seq)

# Output all translations
for i, seq in enumerate(proteins, start=1):
    print(f">protein_{i}\n{seq}\n")

# Optional: save to FASTA file
with open("proteins_E_coli_k1.faa", "w") as f:
    for i, seq in enumerate(proteins, start=1):
        f.write(f">protein_{i}\n{seq}\n")


# Path to your GenBank file
gbk_file_2 = "Fusobacterium_nucleatum_subsp_polymorphum_ATCC_10953.gbk"

# List to store all protein sequences
proteins = []

# Parse GenBank file
for record in SeqIO.parse(gbk_file_2, "genbank"):
    for feature in record.features:
        if feature.type == "CDS" and "translation" in feature.qualifiers:
            protein_seq = feature.qualifiers["translation"][0]
            proteins.append(protein_seq)

# Output all translations
for i, seq in enumerate(proteins, start=1):
    print(f">protein_{i}\n{seq}\n")

# Optional: save to FASTA file
with open("proteins_F_nucleatum_10953.faa", "w") as f:
    for i, seq in enumerate(proteins, start=1):
        f.write(f">protein_{i}\n{seq}\n")