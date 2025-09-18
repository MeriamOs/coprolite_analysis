import pandas as pd

# Load the .bim file
#bim_file = "angsd_bam_trimmed_SE_geno08_filter.bim"
bim_file = "angsd_bam_trimmed_SE_geno08_pre_filter.bim"

df = pd.read_csv(bim_file, sep="\t", header=None, names=["chr", "snp", "cm", "pos", "a1", "a2"], dtype=str)

# Normalize alleles to uppercase
df["a1"] = df["a1"].str.upper()
df["a2"] = df["a2"].str.upper()

# Define damage
is_damage = (
    ((df["a1"] == "C") & (df["a2"] == "T")) |
    ((df["a1"] == "T") & (df["a2"] == "C"))
)

# Extract list of damage SNPs
damage = df.loc[is_damage, "snp"]

# Save to file
damage.to_csv("damage_SE_pre.txt", index=False, header=False)
#damage.to_csv("damage_SE.txt", index=False, header=False)

print(f"Total damage: {len(damage)}")
