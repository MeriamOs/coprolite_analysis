#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=pigs_pca
#SBATCH --time=04:00:00
#SBATCH --mem=8GB
#SBATCH --output PCA.%j.out # CHANGE map1 part each run
#SBATCH --error PCA.%j.err # CHANGE map1 part each run

#prepare input files

#python needs to be version 2
module load Python/2.7.18-gimkl-2020a

vcf_prefix='eager'
#OG8511

#python vcf2eigenstrat.py -v ${vcf_prefix}_samples.vcf -o ${vcf_prefix}


#merge vcf file from Wu and new samples
module load BCFtools/1.16-GCC-11.3.0

#bgzip Wu_modern.vcf
bgzip ${vcf_prefix}_samples.vcf

#bcftools index Wu_modern.vcf.gz
bcftools index ${vcf_prefix}_samples.vcf.gz

#use Wu_modern_new.vcf.gz for PCA calculations
bcftools merge Lapita_merged.vcf.gz ${vcf_prefix}_samples.vcf.gz --missing-to-ref > ${vcf_prefix}_merged.vcf

module purge
module load Python/2.7.18-gimkl-2020a

python vcf2eigenstrat.py -v ${vcf_prefix}_merged.vcf -o ${vcf_prefix}_merged

#run smartpca
module load EIGENSOFT/7.2.1-gimkl-2018b

smartpca -p param_file.txt


#prepare evec file for Rstudio

#new header
header='Individual PC1 PC2 PC3 PC4 Population'
echo $header > header.txt
tr -s " " "\t" < header.txt > ${vcf_prefix}_rstudio.txt

#column selection
awk '{print}' ${vcf_prefix}_merged.evec > ${vcf_prefix}_out.txt
tr -s " " "\t" < ${vcf_prefix}_out.txt > ${vcf_prefix}_out_edit.txt
cat ${vcf_prefix}_out_edit.txt | awk -F "\t" '{print$2, $3, $4, $5, $6, $13}' | tail -n+2 >> ${vcf_prefix}_rstudio.txt

#turn txt into csv file
cat ${vcf_prefix}_rstudio.txt | tr -s '[:blank:]' ',' > ${vcf_prefix}_rstudio.csv
