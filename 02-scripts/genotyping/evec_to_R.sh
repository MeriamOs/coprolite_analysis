vcf_prefix='palaeofaeces_souilmi_merged'

#prepare evec file for Rstudio

#new header
header='Individual PC1 PC2 PC3 PC4 Population'
echo $header > header.txt
tr -s " " "\t" < header.txt > ${vcf_prefix}_rstudio.txt

#column selection
awk '{print}' ${vcf_prefix}.evec > ${vcf_prefix}_out.txt
tr -s " " "\t" < ${vcf_prefix}_out.txt > ${vcf_prefix}_out_edit.txt
cat ${vcf_prefix}_out_edit.txt | awk -F "\t" '{print$2, $3, $4, $5, $6, $13}' | tail -n+2 >> ${vcf_prefix}_rstudio.txt

#turn txt into csv file
cat ${vcf_prefix}_rstudio.txt | tr -s '[:blank:]' ',' > ${vcf_prefix}_rstudio.csv
