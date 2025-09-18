###############################################
### FIXED MUTATIONS PER ARCHAEOLOGICAL SITE ###
###############################################

module load PLINK/1.09b6.16

plink --bfile angsd_bam_trimmed_SE_geno08_no_damage_shared_sites --recode vcf \
    --out angsd_bam_trimmed_SE_geno08_no_damage_shared_sites --dog

python scan_for_fixed_per_site.py

tail -n +2 angsd_bam_trimmed_SE_geno08_no_damage_fixed_snps_by_site.csv \
    | awk -F, 'BEGIN{OFS="\t"} {print $2, $3-1, $3, $1, $4, $5, $6, $7, $8, $9}' \
    > angsd_bam_trimmed_SE_geno08_no_damage_fixed_snps_by_site.bed

module purge
module load BEDTools/2.31.1-GCC-12.3.0

bedtools intersect -a angsd_bam_trimmed_SE_geno08_no_damage_fixed_snps_by_site.bed \
    -b ../pileupcaller_analysis/genes_only_canFam3.1.bed -wa -wb \
    > angsd_bam_trimmed_SE_geno08_no_damage_fixed_snps_by_site_annotated.tsv

### SCRIPT: scan_for_fixed_per_site.py
#import pandas as pd
#
## Load the VCF-like file (assuming no header lines starting with ##)
#vcf = pd.read_csv("angsd_bam_trimmed_SE_geno08_no_damage_shared_sites.vcf", sep="\t")
#
## Extract sample columns (after FORMAT column)
#sample_cols = vcf.columns[9:]
#
## Group sample columns by site name
#site_groups = {}
#for col in sample_cols:
#    site = "_".join(col.split("_")[:-1])  # e.g. Whenua_Hou_post
#    site_groups.setdefault(site, []).append(col)
#
## Function to check if genotypes are fixed for ALT allele
#def get_fixed_alt(genotypes):
#    non_missing = [gt for gt in genotypes if gt != "./."]
#    if len(non_missing) == 0:
#        return None
#    if all(gt == "1/1" for gt in non_missing):
#        return "ALT"
#    return None  # Ignore if fixed for REF or mixed
#
#print(vcf['CHROM'].unique())
#
## Collect results
#output_rows = []
#
#for idx, row in vcf.iterrows():
#    snp_id = f"{row['CHROM']}_{row['POS']}"
#    chrom = row["CHROM"]
#    pos = row["POS"]
#    ref = row["REF"]
#    alt = row["ALT"]
#
#    for site, samples in site_groups.items():
#        genotypes = [row[sample].split(":")[0] for sample in samples]
#        fixed = get_fixed_alt(genotypes)
#        if fixed == "ALT":
#            allele = alt
#            
#            # Count non-missing data for the current site
#            non_missing_count = sum(1 for gt in genotypes if gt != "./.")
#            
#            # Count non-missing data from other sites
#            non_missing_other_sites = 0
#            non_missing_alt_allele = 0
#
#            for other_site, other_samples in site_groups.items():
#                if other_site != site:
#                    for sample in other_samples:
#                        other_genotype = row[sample].split(":")[0]
#                        if other_genotype != "./.":
#                            non_missing_other_sites += 1
#                            if other_genotype == "1/1":
#                                non_missing_alt_allele += 1
#            
#            output_rows.append({
#                "SNP_ID": snp_id,
#                "CHROM": chrom,
#                "POS": pos,
#                "SITE": site,
#                "FIXED_FOR": fixed,
#                "ALLELE": allele,
#                "NON_MISSING_COUNT": non_missing_count,
#                "NON_MISSING_OTHER_SITES": non_missing_other_sites,
#                "NON_MISSING_ALT_ALLELE": non_missing_alt_allele
#            })
#
## Convert to DataFrame and save
#results_df = pd.DataFrame(output_rows)
#results_df.to_csv("angsd_bam_trimmed_SE_geno08_no_damage_fixed_snps_by_site.csv", index=False)
#
#print("✅ Results saved to 'angsd_bam_trimmed_SE_geno08_no_damage_fixed_snps_by_site.csv'")
