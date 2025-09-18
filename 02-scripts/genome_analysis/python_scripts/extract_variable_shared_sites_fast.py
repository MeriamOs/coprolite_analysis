from collections import Counter

def filter_shared_variable_sites(
    geno_file="angsd_bam_trimmed_SE_geno08_no_damage.geno",
    snp_file="angsd_bam_trimmed_SE_geno08_no_damage.snp",
    output_geno="angsd_bam_trimmed_SE_geno08_no_damage_shared_sites.geno",
    output_snp="angsd_bam_trimmed_SE_geno08_no_damage_shared_sites.snp"
):
    variable_indices = []
    variable_lines = []

    # Step 1: Identify shared variable SNPs
    with open(geno_file) as f:
        for idx, line in enumerate(f):
            line = line.strip()
            alleles = [g for g in line if g != '9']
            if len(set(alleles)) > 1:
                counts = Counter(alleles)
                sorted_counts = sorted(counts.values(), reverse=True)
                if len(sorted_counts) >= 2 and sorted_counts[1] >= 2:
                    variable_indices.append(idx)
                    variable_lines.append(line)

    # Step 2: Save filtered .geno
    with open(output_geno, 'w') as g_out:
        g_out.write('\n'.join(variable_lines) + '\n')

    # Step 3: Save corresponding .snp lines using a set for fast lookup
    variable_set = set(variable_indices)
    with open(snp_file) as s_in, open(output_snp, 'w') as s_out:
        for idx, line in enumerate(s_in):
            if idx in variable_set:
                s_out.write(line)

    print(f"Saved {len(variable_lines)} shared variable SNPs.")
    print(f"  - GENO: {output_geno}")
    print(f"  - SNP:  {output_snp}")

# Run the function
filter_shared_variable_sites()
