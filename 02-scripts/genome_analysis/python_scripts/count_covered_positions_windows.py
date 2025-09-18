import pandas as pd

# Parameters
vcf_file = "angsd_bam_trimmed_SE_no_post_04_haploid.vcf" 
window_size = 10000

# Read VCF file (skip header lines)
vcf_data = pd.read_csv(vcf_file, sep='\t', comment='#', header=None,
                       usecols=[0, 1], names=['CHROM', 'POS'])

# Group by CHROM
window_counts = []

for chrom, group in vcf_data.groupby('CHROM'):
    # Sort positions just in case
    positions = group['POS'].sort_values().values
    max_pos = positions[-1]

    for start in range(0, max_pos + 1, window_size):
        end = start + window_size
        count = ((positions >= start) & (positions < end)).sum()
        window_counts.append({
            'CHROM': chrom,
            'START': start,
            'END': end,
            'N_SITES': count
        })

# Convert to DataFrame and save
result_df = pd.DataFrame(window_counts)
result_df.to_csv('positions_counts_per_10kb_window.tsv', sep='\t', index=False)
