import pandas as pd
import matplotlib.pyplot as plt

# Load the full file
df = pd.read_csv("angsd_bam_trimmed_SE_all_05_haploid_custom.frq", delim_whitespace=True)
#df = pd.read_csv("angsd_bam_trimmed_SE_no_post_04_haploid_custom.frq", delim_whitespace=True)

print(f"Loaded {len(df)} rows from the .frq file")

# Extract alt allele frequency from ALLELE2:FREQ
alt_freq = df['ALLELE2:FREQ'].str.split(':').str[1].astype(float)
n_chr = df['N_CHR'].astype(int)

# Calculate alternative allele count
alt_count = (alt_freq * n_chr).round().astype(int)

# Count SFS (Site Frequency Spectrum), including fixed reference (0) and fixed alt (n_chr)
sfs = alt_count.value_counts().sort_index()

print("\nAlt Allele Count\tNumber of Sites")
for count, num_snps in sfs.items():
    print(f"{count}\t{num_snps}")

# Plot
plt.bar(sfs.index, sfs.values, color="purple")
plt.xlabel("Alternative Allele Count")
plt.ylabel("Number of Sites")
plt.title("Site Frequency Spectrum (including fixed sites)")
plt.tight_layout()

plt.savefig("angsd_bam_trimmed_SE_all_05_haploid_alt_sfs_plot.png", dpi=300)
#plt.savefig("angsd_bam_trimmed_SE_no_post_04_haploid_alt_sfs_plot.png", dpi=300)
