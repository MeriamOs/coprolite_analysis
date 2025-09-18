import csv

#input_vcf = "../angsd_bam_trimmed_SE_no_post_04_haploid.vcf"
#input_vcf = "angsd_bam_trimmed_SE_no_post_04_haploid_LB.vcf"
#input_vcf = "angsd_bam_trimmed_SE_no_post_04_haploid_WH.vcf"
input_vcf = "angsd_bam_trimmed_SE_no_post_04_haploid_KH.vcf"

#output_frq = "angsd_bam_trimmed_SE_no_post_04_haploid_custom.frq"
#output_frq = "angsd_bam_trimmed_SE_no_post_04_haploid_LB_custom.frq"
#output_frq = "angsd_bam_trimmed_SE_no_post_04_haploid_WH_custom.frq"
output_frq = "angsd_bam_trimmed_SE_no_post_04_haploid_KH_custom.frq"

with open(input_vcf, "r") as fin, open(output_frq, "w", newline="") as fout:
    reader = csv.reader(fin, delimiter="\t")
    writer = csv.writer(fout, delimiter="\t")

    header_written = False

    for row in reader:
        if row[0].startswith("##"):
            continue
        elif row[0] == "#CHROM":
            sample_names = row[9:]
            if not header_written:
                writer.writerow(["CHROM", "POS", "N_ALLELES", "N_CHR", "ALLELE1:FREQ", "ALLELE2:FREQ"])
                header_written = True
            continue

        chrom, pos, _, ref, alt = row[0], row[1], row[2], row[3], row[4]
        genotypes = row[9:]

        ref_count = sum(1 for g in genotypes if g == "0")
        alt_count = sum(1 for g in genotypes if g == "1")
        total = ref_count + alt_count

        if total == 0:
            continue

        ref_freq = ref_count / total
        alt_freq = alt_count / total

        writer.writerow([
            chrom,
            pos,
            "1",
            total,
            f"{ref}:{ref_freq:.3f}",  # REF is ALLELE1
            f"{alt}:{alt_freq:.3f}"   # ALT is ALLELE2
        ])
