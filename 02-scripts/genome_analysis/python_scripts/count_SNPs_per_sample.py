import pysam
from collections import defaultdict

# Input VCF file
vcf_path = "../angsd_bam_trimmed_SE_no_post_04_haploid.vcf"

# Open VCF
vcf = pysam.VariantFile(vcf_path)

# Initialize counters per sample
sample_stats = {sample: {"covered": 0, "alt": 0} for sample in vcf.header.samples}

# Iterate through records
for rec in vcf.fetch():
    for sample in vcf.header.samples:
        sample_data = rec.samples[sample]
        gt = sample_data.get("GT")

        # Skip if no genotype
        if gt is None or all(allele is None for allele in gt):
            continue

        # Covered site
        sample_stats[sample]["covered"] += 1

        # ALT site if genotype includes any non-zero allele
        if any(allele not in (None, 0) for allele in gt):
            sample_stats[sample]["alt"] += 1

# Print results
print("Sample\tCovered_Sites\tAlt_Allele_Sites")
for sample, stats in sample_stats.items():
    print(f"{sample}\t{stats['covered']}\t{stats['alt']}")

with open("sample_counts.tsv", "w") as fout:
    fout.write("Sample\tCovered_Sites\tAlt_Allele_Sites\n")
    for sample, stats in sample_stats.items():
        fout.write(f"{sample}\t{stats['covered']}\t{stats['alt']}\n")