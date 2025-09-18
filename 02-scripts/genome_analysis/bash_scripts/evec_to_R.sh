# Script to prepare evec file for Rstudio 

# NOTE: change vcf prefix below
vcf_prefix='angsd_bam_trimmed_geno02_transversions_only'

# new header
header='Individual PC1 PC2 PC3 PC4 Population'
echo $header > header.txt
tr -s " " "\t" < header.txt > ${vcf_prefix}_rstudio.txt

# column selection
awk '{print}' ${vcf_prefix}.evec > ${vcf_prefix}_out.txt
tr -s " " "\t" < ${vcf_prefix}_out.txt > ${vcf_prefix}_out_edit.txt
sed 's/^\t//' ${vcf_prefix}_out_edit.txt > ${vcf_prefix}_out_edit2.txt
cat ${vcf_prefix}_out_edit2.txt | awk -F "\t" '{print$1, $2, $3, $4, $5, $12}' | tail -n+2 >> ${vcf_prefix}_rstudio.txt

# turn txt into csv file
cat ${vcf_prefix}_rstudio.txt | tr -s '[:blank:]' ',' > ${vcf_prefix}_rstudio.csv
rm ${vcf_prefix}_out.txt ${vcf_prefix}_out_edit.txt ${vcf_prefix}_out_edit2.txt