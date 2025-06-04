#downloaded Souilmi plink files from  https://doi.org/10.25909/25885747, and renamed them to shorten file names
#souilmi2024_bergstrom2022_bergstrom2020_zhang2020_plassais2019_finalset.bergsgrom2022-analyses-loci.tv.c2
#to souilmi2024.*

module load PLINK/1.09b6.16

#checked for pseudohaploid data
plink --bfile original/souilmi2024_bergstrom2022_bergstrom2020_zhang2020_plassais2019_finalset.bergsgrom2022-analyses-loci.tv.c2 \
  --dog --het --out heterozygosity

#Pull out ancient samples, except dingoes to remove
awk '$6 == 1 && $2 !~ /_ph$/ {print $1, $2}' heterozygosity.het > pseudohaploid_samples.txt

#Pull out *unkown* samples, and ancient dingoes with *_dp 
awk '$2 ~ /_dp/ || $1 ~ /nknown/ {print $1, $2}' heterozygosity.het > dp_unknown_samples.txt

cat pseudohaploid_samples.txt dp_unknown_samples.txt > souilmi2024_to_remove.txt
#also added fox and dhole, to minimize outgroup (gene flow between groups)

plink --bfile souilmi2024 --remove souilmi2024_to_remove.txt --make-bed \
    --out souilmi2024_filtered --dog
#Changed the population names in *.fam, e.g. dingo, dingo_ancient, and wolfLocation

cd ..
cat souilmi_files/souilmi2024_filtered > merge_list.txt

cp ../eager_CanFam3.1/genotyping/pileupcaller.single.* .

bash snp_to_plink.sh 
#snp_to_plink.sh
#!/bin/bash
#module load EIGENSOFT/7.2.1-gimkl-2018b
#in1='pileupcaller.single' 
#
#convertf -p <(echo "genotypename:	${in1}.geno
#snpname:	${in1}.snp
#indivname:	${in1}.ind
#outputformat:	PACKEDPED
#genotypeoutname:	${in1}.bed
#snpoutname:	${in1}.bim
#indivoutname:	${in1}.fam")

#change *.fam file to
#Kuri MS11669 0 0 2 1
#Kuri MS11670 0 0 1 1
#Kuri MS11673 0 0 2 1
#Kuri MS11674 0 0 2 1
#Kuri MS11675 0 0 1 1
#Kuri MS11676 0 0 2 1
#Kuri MS11677 0 0 2 1
#Kuri MS11678 0 0 2 1
#Kuri MS11679 0 0 1 1
#Kuri MS11683 0 0 1 1
#Kuri MS11684 0 0 1 1
#Kuri MS11686 0 0 2 1
#Kuri MS11770 0 0 1 1
#Kuri MS11771 0 0 1 1
#Kuri MS11774 0 0 1 1
#Kuri MS11775 0 0 1 1

plink --bfile pileupcaller.single --merge-list merge_list.txt \
    --make-bed --out merged_kuri_souilmi --dog --allow-no-sex
#489 samples

#### CREATE SUBSET WITHOUT OTHER CANINE SPECIES ####
awk '{print $1, $2}' merged_kuri_souilmi.fam > souilmi_kuri_nowolves_to_keep.txt 
#removed all wolves, coyotes, jackals, and african hunting dog

plink --bfile merged_kuri_souilmi --keep souilmi_kuri_nowolves_to_keep.txt \
    --make-bed --out merged_kuri_souilmi_nowolves --dog
#409 samples

#### CREATE SUBSET WITH EAST ASIAN DESCENDENT DOGS ####
awk '{print $1, $2}' souilmi_kuri_nowolves_to_keep.txt > souilmi_kuri_subset_to_keep.txt
#kept only dingoes, NGSD, East Asian village and indigenous dogs, Xiasi, Jindo, and Chow chows

plink --bfile merged_kuri_souilmi --keep souilmi_kuri_subset_to_keep.txt \
    --make-bed --out merged_kuri_souilmi_subset --dog
#92 samples

module load EIGENSOFT/7.2.1-gimkl-2018b

#create list to ignore kuri and ancient dingoes for PC computation 
awk '{print $1}' merged_kuri_souilmi.fam | grep -v -E '^(Kuri|Dingo_ancient)$' | sort | uniq > population_list_global.txt
awk '{$NF = $1} 1' merged_kuri_souilmi.fam > merged_kuri_souilmi_modified.fam

awk '{print $1}' merged_kuri_souilmi_subset.fam | grep -v -E '^(Kuri|Dingo_ancient)$' | sort | uniq > population_list_subset.txt
awk '{$NF = $1} 1' merged_kuri_souilmi_subset.fam > merged_kuri_souilmi_subset_modified.fam

awk '{print $1}' merged_kuri_souilmi_nowolves.fam | grep -v -E '^(Kuri|Dingo_ancient)$' | sort | uniq > population_list_no_wolves.txt
awk '{$NF = $1} 1' merged_kuri_souilmi_nowolves.fam > merged_kuri_souilmi_nowolves_modified.fam

######################
#### RUN SMARTPCA ####
######################

bash smartpca.sh

#smartpca.sh
#module load EIGENSOFT/7.2.1-gimkl-2018b
#smartpca -p smartpca_global.param
#smartpca -p smartpca_nowolves.param
#smartpca -p smartpca_subset.param

#smartpca_subset.param
#genotypename:  merged_kuri_souilmi_subset.bed
#snpname:   merged_kuri_souilmi_subset.bim
#indivname:   merged_kuri_souilmi_subset_modified.fam
#evecoutname:  merged_kuri_souilmi_subset.evec
#evaloutname:   merged_kuri_souilmi_subset.eval
#poplistname: population_list_subset.txt
#lsqproject: YES
#newshrink: NO
#missingmode: NO
#inbreed: YES
#numchrom: 38

#smartpca_global.param
#genotypename:  merged_kuri_souilmi.bed
#snpname:   merged_kuri_souilmi.bim
#indivname:   merged_kuri_souilmi_modified.fam
#evecoutname:  merged_kuri_souilmi.evec
#evaloutname:   merged_kuri_souilmi.eval
#poplistname: population_list_global.txt
#lsqproject: YES
#newshrink: NO
#missingmode: NO
#inbreed: YES
#numchrom: 38
#nthreads: 32

#smartpca_nowolves.param
#genotypename:  merged_kuri_souilmi_nowolves.bed
#snpname:   merged_kuri_souilmi_nowolves.bim
#indivname:   merged_kuri_souilmi_nowolves_modified.fam
#evecoutname:  merged_kuri_souilmi_nowolves.evec
#evaloutname:   merged_kuri_souilmi_nowolves.eval
#poplistname: population_list_no_wolves.txt
#lsqproject: YES
#newshrink: NO
#missingmode: NO
#inbreed: YES
#numchrom: 38

bash evec_to_R.sh

### evec_to_R.sh
#vcf_prefix='merged_kuri_souilmi_subset'
#
##new header
#header='Individual PC1 PC2 PC3 PC4 Population'
#echo $header > header.txt
#tr -s " " "\t" < header.txt > ${vcf_prefix}_rstudio.txt
#
##column selection
#awk '{print}' ${vcf_prefix}.evec > ${vcf_prefix}_out.txt
#tr -s " " "\t" < ${vcf_prefix}_out.txt > ${vcf_prefix}_out_edit.txt
#sed 's/^\t//' ${vcf_prefix}_out_edit.txt > ${vcf_prefix}_out_edit2.txt
#cat ${vcf_prefix}_out_edit2.txt | awk -F "\t" '{print$1, $2, $3, $4, $5, $12}' | tail -n+2 >> ${vcf_prefix}_rstudio.txt
#
##turn txt into csv file
#cat ${vcf_prefix}_rstudio.txt | tr -s '[:blank:]' ',' > ${vcf_prefix}_rstudio.csv
#rm ${vcf_prefix}_out.txt ${vcf_prefix}_out_edit.txt ${vcf_prefix}_out_edit2.txt ${vcf_prefix}_rstudio.txt