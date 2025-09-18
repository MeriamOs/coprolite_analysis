module load PLINK/1.09b6.16

#####################
#### GENOME WIDE ####
#####################

#includes autosomes and X chromosome
#use dataset that is not filtered with maf to retain non variable positions
plink --bfile ../angsd_bam_trimmed_SE_geno08_pre --biallelic-only \
       --dog --make-bed --out angsd_bam_trimmed_SE_no_post --snps-only

plink --bfile angsd_bam_trimmed_SE_no_post --dog \
       --recode vcf --out angsd_bam_trimmed_SE_no_post

module load BCFtools/1.19-GCC-11.3.0

nano samples_bcftools.txt

### FILE: samples_bcftools.txt
#Whenua_Hou_post_MS11669 F
#Whenua_Hou_post_MS11670 M
#Whenua_Hou_post_MS11673 F 
#Whenua_Hou_pre_MS11674 F 
#Whenua_Hou_pre_MS11675 M 
#Whenua_Hou_pre_MS11676 F 
#Whenua_Hou_pre_MS11677 F 
#Whenua_Hou_pre_MS11678 F 
#Long_Bay_MS11679 M 
#Long_Bay_MS11683 M 
#Long_Bay_MS11684 M 
#Long_Bay_MS11686 F 
#Kahukura_MS11770 M 
#Kahukura_MS11771 M 
#Kahukura_MS11774 M 
#Kahukura_MS11775 M 

bcftools +fixploidy angsd_bam_trimmed_SE_no_post.vcf -- -s samples_bcftools.txt -f 1 \
       | bcftools view -i 'F_MISSING <= 0.4' \
       | bcftools view -i 'REF!="C" || ALT!="T"' \
       | bcftools view -i 'REF!="T" || ALT!="C"' -Ov -o angsd_bam_trimmed_SE_no_post_04_haploid.vcf
#138,738,633 variants remaining

vcftools --vcf angsd_bam_trimmed_SE_no_post_04_haploid.vcf --window-pi 10000 \
       --haploid --out angsd_bam_trimmed_SE_no_post_04_haploid

vcftools --vcf angsd_bam_trimmed_SE_no_post_04_haploid.vcf --TajimaD 10000 \
       --haploid --out angsd_bam_trimmed_SE_no_post_04_haploid

module load Python/3.11.6-foss-2023a

python count_covered_positions_windows.py

### SCRIPT: count_covered_positions_windows.py
#
#import pandas as pd
#
## Parameters
#vcf_file = "angsd_bam_trimmed_SE_no_post_04_haploid.vcf" 
#window_size = 10000
#
## Read VCF file (skip header lines)
#vcf_data = pd.read_csv(vcf_file, sep='\t', comment='#', header=None,
#                       usecols=[0, 1], names=['CHROM', 'POS'])
#
## Group by CHROM
#window_counts = []
#
#for chrom, group in vcf_data.groupby('CHROM'):
#    # Sort positions just in case
#    positions = group['POS'].sort_values().values
#    max_pos = positions[-1]
#
#    for start in range(0, max_pos + 1, window_size):
#        end = start + window_size
#        count = ((positions >= start) & (positions < end)).sum()
#        window_counts.append({
#            'CHROM': chrom,
#            'START': start,
#            'END': end,
#            'N_SITES': count
#        })
#
## Convert to DataFrame and save
#result_df = pd.DataFrame(window_counts)
#result_df.to_csv('positions_counts_per_10kb_window.tsv', sep='\t', index=False)


###########################
#### WINDOW ANNOTATION ####
###########################

#downloaded annotation *.gff from GenBank, GCF_000002285.3

#cp genomic.gff genome_canFam3.1_annotations.gff

awk '$3 == "gene"' genome_canFam3.1_annotations.gff > genes_only_canFam3.1_annotations.gff3
awk 'FNR==NR {map[$3]=$1; next} {if($1 in map) $1=map[$1]; print}' OFS="\t" \
    chromosome_sequence_mapping.tsv genes_only_canFam3.1_annotations.gff3 \
    > genes_only_canFam3.1_annotations_names_fixed.gff3

awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, $9}' \
        genes_only_canFam3.1_annotations_names_fixed.gff3 \
        > genes_only_canFam3.1.bed

module load BEDTools/2.31.1-GCC-12.3.0

awk 'NR>1 {print $1"\t"($2)"\t"$3"\t"$4"\t"$5}' \
    angsd_bam_trimmed_SE_no_post_04_haploid.windowed.pi \
    > angsd_bam_trimmed_SE_no_post_04_haploid.windowed.bed

awk 'NR>1 {print $1"\t"($2)"\t"($2+9999)"\t"$3"\t"$4}' \
    angsd_bam_trimmed_SE_no_post_04_haploid.Tajima.D \
    > angsd_bam_trimmed_SE_no_post_04_haploid.Tajima.D.bed

bedtools intersect -a angsd_bam_trimmed_SE_no_post_04_haploid.windowed.bed \
    -b ../../genotyping/pileupcaller_analysis/genes_only_canFam3.1.bed \
    -wa -wb > angsd_bam_trimmed_SE_no_post_04_haploid_annotated.bed

bedtools intersect -a angsd_bam_trimmed_SE_no_post_04_haploid.Tajima.D.bed \
    -b ../../genotyping/pileupcaller_analysis/genes_only_canFam3.1.bed \
    -wa -wb > angsd_bam_trimmed_SE_no_post_04_haploid_annotated.Tajima.D.bed

#########################
#### MITOGENOMES VCF ####
#########################

#in eager/genotyping folder
bcftools merge -Oz -o kuri_all_mtdna.vcf.gz MS11669.unifiedgenotyper.vcf.gz MS11670.unifiedgenotyper.vcf.gz \
       MS11673.unifiedgenotyper.vcf.gz MS11674.unifiedgenotyper.vcf.gz MS11675.unifiedgenotyper.vcf.gz \
       MS11676.unifiedgenotyper.vcf.gz MS11677.unifiedgenotyper.vcf.gz MS11678.unifiedgenotyper.vcf.gz \
       MS11679.unifiedgenotyper.vcf.gz MS11683.unifiedgenotyper.vcf.gz MS11684.unifiedgenotyper.vcf.gz \
       MS11686.unifiedgenotyper.vcf.gz MS11770.unifiedgenotyper.vcf.gz MS11771.unifiedgenotyper.vcf.gz \
       MS11774.unifiedgenotyper.vcf.gz MS11775.unifiedgenotyper.vcf.gz

gunzip kuri_all_mtdna.vcf.gz

#mask variable number of tandem repeats (VNTR) region - 16130-16430
bcftools +fixploidy kuri_all_mtdna.vcf -- -s samples_kuri_bcftools.txt -f 1 \
       | bcftools view -e 'POS>=16130 && POS<=16430' \
       | bcftools view -Ov -o kuri_all_mtdna_haploid.vcf

vcftools --vcf kuri_all_mtdna_haploid.vcf --window-pi 1000 \
       --haploid --out kuri_all_mtdna_pi

awk 'NR > 1 {sum += $5} END {print "Mean PI:", sum / (NR - 1)}' \
       kuri_all_mtdna_pi.windowed.pi 
