module load PLINK/1.09b6.16

#####################
#### Genome wide ####
#####################

#includes autosomes and X chromosome

plink --bfile ../angsd_bam_trimmed_SE_geno08_pre --biallelic-only \
       --dog --make-bed --out angsd_bam_trimmed_SE_no_post --snps-only

plink --bfile angsd_bam_trimmed_SE_no_post --dog \
       --recode vcf --out angsd_bam_trimmed_SE_no_post

module load BCFtools/1.19-GCC-11.3.0

bcftools +fixploidy angsd_bam_trimmed_SE_no_post.vcf -- -s samples_bcftools.txt -f 1 \
       | bcftools view -i 'F_MISSING <= 0.4' \
       | bcftools view -i 'REF!="C" || ALT!="T"' \
       | bcftools view -i 'REF!="T" || ALT!="C"' -Ov -o angsd_bam_trimmed_SE_no_post_04_haploid.vcf
#138,738,633 variants remaining

vcftools --vcf angsd_bam_trimmed_SE_no_post_04_haploid.vcf --window-pi 10000 \
       --haploid --out angsd_bam_trimmed_SE_no_post_04_haploid

vcftools --vcf angsd_bam_trimmed_SE_no_post_04_haploid.vcf --TajimaD 10000 \
       --haploid --out angsd_bam_trimmed_SE_no_post_04_haploid

###########################
#### window annotation ####
###########################

#downloaded annotation *.gff from GenBank, GCF_000002285.3

#cp genomic.gff genome_canFam3.1_annotations.gff

#awk '$3 == "gene"' genome_canFam3.1_annotations.gff > genes_only_canFam3.1_annotations.gff3
#awk 'FNR==NR {map[$3]=$1; next} {if($1 in map) $1=map[$1]; print}' OFS="\t" \
#    chromosome_sequence_mapping.tsv genes_only_canFam3.1_annotations.gff3 \
#    > genes_only_canFam3.1_annotations_names_fixed.gff3

#awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, $9}' \
#        genes_only_canFam3.1_annotations_names_fixed.gff3 \
#        > genes_only_canFam3.1.bed

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
#### mitogenomes vcf ####
#########################

#mask variable number of tandem repeats (VNTR) region - 16130-16430
bcftools +fixploidy kuri_all_mtdna.vcf -- -s samples_kuri_bcftools.txt -f 1 \
       | bcftools view -e 'POS>=16130 && POS<=16430' \
       | bcftools view -Ov -o kuri_all_mtdna_haploid.vcf

vcftools --vcf kuri_all_mtdna_haploid.vcf --window-pi 1000 \
       --haploid --out kuri_all_mtdna_pi

awk 'NR > 1 {sum += $5} END {print "Mean PI:", sum / (NR - 1)}' \
       kuri_all_mtdna_pi.windowed.pi 
