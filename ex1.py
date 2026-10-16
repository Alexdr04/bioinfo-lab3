from collections import Counter
import matplotlib.pyplot as plt

FASTA_FILE = "GCF_002605485.1_ASM260548v1_cds_from_genomic.fna"

codon_table = {
    "UUU": "Phenylalanine", "UUC": "Phenylalanine", "UUA": "Leucine", "UUG": "Leucine",
    "CUU": "Leucine", "CUC": "Leucine", "CUA": "Leucine", "CUG": "Leucine",
    "AUU": "Isoleucine", "AUC": "Isoleucine", "AUA": "Isoleucine", "AUG": "Methionine",
    "GUU": "Valine", "GUC": "Valine", "GUA": "Valine", "GUG": "Valine",
    "UCU": "Serine", "UCC": "Serine", "UCA": "Serine", "UCG": "Serine",
    "CCU": "Proline", "CCC": "Proline", "CCA": "Proline", "CCG": "Proline",
    "ACU": "Threonine", "ACC": "Threonine", "ACA": "Threonine", "ACG": "Threonine",
    "GCU": "Alanine", "GCC": "Alanine", "GCA": "Alanine", "GCG": "Alanine",
    "UAU": "Tyrosine", "UAC": "Tyrosine", "UAA": "Stop", "UAG": "Stop",
    "CAU": "Histidine", "CAC": "Histidine", "CAA": "Glutamine", "CAG": "Glutamine",
    "AAU": "Asparagine", "AAC": "Asparagine", "AAA": "Lysine", "AAG": "Lysine",
    "GAU": "Aspartic acid", "GAC": "Aspartic acid", "GAA": "Glutamic acid", "GAG": "Glutamic acid",
    "UGU": "Cysteine", "UGC": "Cysteine", "UGA": "Stop", "UGG": "Tryptophan",
    "CGU": "Arginine", "CGC": "Arginine", "CGA": "Arginine", "CGG": "Arginine",
    "AGU": "Serine", "AGC": "Serine", "AGA": "Arginine", "AGG": "Arginine",
    "GGU": "Glycine", "GGC": "Glycine", "GGA": "Glycine", "GGG": "Glycine",
}

sequence = ""
with open(FASTA_FILE) as f:
    for line in f:
        if not line.startswith(">"):
            sequence += line.strip().upper()

counts = Counter()
for i in range(len(sequence) - 2):
    codon = sequence[i:i + 3]
    if set(codon) <= set("ACGT"):
        counts[codon] += 1

total = sum(counts.values())
top10 = counts.most_common(10)

print(f"Sequence length: {len(sequence)} nt, codons counted: {total}\n")
print(f"{'#':<3}{'DNA':<6}{'RNA':<6}{'Count':<8}{'Freq %':<8}Amino acid")
for rank, (codon, count) in enumerate(top10, start=1):
    rna = codon.replace("T", "U")
    aa = codon_table[rna]
    print(f"{rank:<3}{codon:<6}{rna:<6}{count:<8}{count / total * 100:<8.2f}{aa}")

labels = [f"{c}\n{codon_table[c.replace('T', 'U')]}" for c, _ in top10]
plt.bar(labels, [n for _, n in top10])
plt.title("Top 10 codons")
plt.ylabel("Count")
plt.tight_layout()
plt.show()
