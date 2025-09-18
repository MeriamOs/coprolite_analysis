#downloaded annotation *.gff from GenBank, GCF_000002285.3

#cp genomic.gff genome_canFam3.1_annotations.gff

#awk '$3 == "gene"' genome_canFam3.1_annotations.gff > genes_only_canFam3.1_annotations.gff3
#awk 'FNR==NR {map[$3]=$1; next} {if($1 in map) $1=map[$1]; print}' OFS="\t" chromosome_sequence_mapping.tsv genes_only_canFam3.1_annotations.gff3 > genes_only_canFam3.1_annotations_names_fixed.gff3
#awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, $9}' genes_only_canFam3.1_annotations_names_fixed.gff3 > genes_only_canFam3.1.bed

###########################################
### PER GENE VARIABLE SHARED SNP COUNTS ###
###########################################

awk '{start = $4 - 1; print $2 "\t" start "\t" $4 "\t" $1}' \
    angsd_bam_trimmed_SE_geno08_no_damage_shared_sites.snp \
    > angsd_bam_trimmed_SE_geno08_no_damage_shared_sites_positions.bed

module load BEDTools/2.31.1-GCC-12.3.0

bedtools intersect -a angsd_bam_trimmed_SE_geno08_no_damage_shared_sites_positions.bed \
    -b ../pileupcaller_analysis/genes_only_canFam3.1.bed -wa -wb \
    > angsd_bam_trimmed_SE_geno08_no_damage_shared_snps_in_genes.tsv

module purge
module load Python/3.11.6-foss-2023a

python count_snps_per_gene.py 
#Warning: Inconsistent gene lengths found for the following genes:
#  DEFB110: lengths found = [30400, 5121] (using most common: 30400)
#  OASL: lengths found = [13014, 13718] (using most common: 13014)
#  TRNAR-UCU: lengths found = [88, 92, 75] (using most common: 75)
#  TRNAY-GUA: lengths found = [91, 89, 92] (using most common: 89)

### SCRIPT: count_snps_per_gene.py
#import pandas as pd
#import re
#from collections import Counter
#
## Read the SNP data
#df = pd.read_csv("angsd_bam_trimmed_SE_geno08_no_damage_shared_snps_in_genes.tsv", sep="\t", header=None)
#
## Extract relevant columns
#df['gene_annotations'] = df[7]
#df['gene_name'] = df['gene_annotations'].apply(lambda x: re.search(r'Name=([^;]+)', x).group(1) if pd.notnull(x) else None)
#df['gene_start'] = df[5]
#df['gene_end'] = df[6]
#df['gene_length'] = df['gene_end'] - df['gene_start'] + 1
#
## Drop rows with missing gene names
#df = df.dropna(subset=['gene_name'])
#
## Check for inconsistent gene lengths
#length_issues = {}
#
## Create a dictionary to store the most common length per gene
#gene_lengths = {}
#
#for gene, group in df.groupby('gene_name'):
#    unique_lengths = group['gene_length'].unique()
#    if len(unique_lengths) > 1:
#        length_issues[gene] = unique_lengths
#        # Use the most frequent length
#        common_length = group['gene_length'].mode().iloc[0]
#    else:
#        common_length = unique_lengths[0]
#    gene_lengths[gene] = common_length
#
## Report any issues
#if length_issues:
#    print("Warning: Inconsistent gene lengths found for the following genes:")
#    for gene, lengths in length_issues.items():
#        print(f"  {gene}: lengths found = {list(lengths)} (using most common: {gene_lengths[gene]})")
#
## Count SNPs per gene
#snp_counts = df['gene_name'].value_counts()
#
## Combine counts and gene lengths into a DataFrame
#result = pd.DataFrame({
#    'gene_name': snp_counts.index,
#    'snp_count': snp_counts.values,
#    'gene_length': snp_counts.index.map(gene_lengths)
#})
#
## Calculate mutation rate
#result['mutation_rate'] = result['snp_count'] / result['gene_length']
#
## Save to file
#result.to_csv("angsd_bam_trimmed_SE_geno08_no_damage_shared_snp_rate_per_gene.csv", index=False)
#
## Print result
#print(result)
