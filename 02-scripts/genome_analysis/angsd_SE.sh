#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=angsd
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --output angsd.%j.out # CHANGE map1 part each run
#SBATCH --error angsd.%j.err # CHANGE map1 part each runmodule purge 

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/effective_population

module purge

module load angsd/0.935-GCC-9.2.0
#module load SAMtools

#############################
### CREATE ANGSD DATASETS ###
#############################

#mkdir bam_files
#cd bam_files
#ln -s ../../eager_CanFam3.1/trimmed_bam/* .
#cd ..

#ls bam_files/* > file_list_angsd_SE.txt

angsd -bam file_list_angsd.txt \
    -GL 1 \
    -doCounts 1 \
    -doMajorMinor 1 \
    -doMaf 1 \
    -minMapQ 30 \
    -minQ 30 \
    -uniqueOnly 1 \
    -doHaploCall 1 \
    -out angsd_bam_trimmed_SE \
    -P 16

haploToPlink angsd_bam_trimmed_SE.haplo.gz angsd_bam_trimmed_SE 
#this outputs .tped/.tfam files, that need to be converted to binary plink files

module load PLINK/1.09b6.16
#.tped file had too many variants 
#(Error: PLINK does not support more than 2^31 - 3 variants = 2,147,483,645 variants), so had to fiter it first

awk '($1 ~ /^([1-9]|[12][0-9]|3[0-8]|X)$/)' angsd_bam_trimmed_SE.tped \
        > angsd_bam_trimmed_SE_awk.tped
cp angsd_bam_trimmed_SE.tfam angsd_bam_trimmed_SE_awk.tfam

#Remove positions that have more than 80% missing data 
#meaning that at least 4/16 samples have data
awk '{ n_miss = 0; total = 0;
  for (i=5; i<=NF; i++) { total++; if ($i == "N") n_miss++;}
  miss_rate = n_miss / total;
  if (miss_rate <= 0.8) print;}' angsd_bam_trimmed_SE_awk.tped \
        > angsd_bam_trimmed_SE_geno08.tped
cp angsd_bam_trimmed_SE.tfam angsd_bam_trimmed_SE_geno08.tfam

plink --tfile angsd_bam_trimmed_SE_geno08 --missing-genotype N \
        --make-bed --out angsd_bam_trimmed_SE_geno08 --dog
#1,698,386,048 variants remaining

#check for no. covered positions per individual
plink --bfile angsd_bam_trimmed_SE_geno08 --missing --dog --out angsd_bam_trimmed_SE_geno08_missingness

plink --bfile angsd_bam_trimmed_SE_geno08 --missing-genotype N \
        --make-bed --dog --snps-only --biallelic-only --maf 0.01 \
        --out angsd_bam_trimmed_SE_geno08_filter  
#51,031,252 variants remaining, these are all variable as per --maf 0.01

module load Python/3.11.6-foss-2023a

## identify and remove damage - I have only removed C-to-T and T-to-C changes considering single stranded protocol
python identify_CtoT.py
#Total identified: 24,151,308

### SCRIPT: identify_CtoT.py
#import pandas as pd
#
## Load the .bim file
#bim_file = "angsd_bam_trimmed_SE_geno08_pre_filter.bim"
#
#df = pd.read_csv(bim_file, sep="\t", header=None, names=["chr", "snp", "cm", "pos", "a1", "a2"], dtype=str)
#
## Normalize alleles to uppercase
#df["a1"] = df["a1"].str.upper()
#df["a2"] = df["a2"].str.upper()
#
## Define damage
#is_damage = (
#    ((df["a1"] == "C") & (df["a2"] == "T")) |
#    ((df["a1"] == "T") & (df["a2"] == "C"))
#)
#
## Extract list of damage SNPs
#damage = df.loc[is_damage, "snp"]
#
## Save to file
#damage.to_csv("damage_SE.txt", index=False, header=False)
#
#print(f"Total damage: {len(damage)}")

plink --bfile angsd_bam_trimmed_SE_geno08_filter --missing-genotype N \
        --exclude damage_SE.txt --make-bed --dog \
        --out angsd_bam_trimmed_SE_geno08_no_damage
#26,879,944 variants remaining

### > now we have a dataset that is quality filtered, only variable sites and has damage removed < ###

## identify and remove transitions
#python identify_transitions.py

#plink --bfile angsd_bam_trimmed_SE_geno08_filter --missing-genotype N \
#        --exclude transitions_SE.txt --make-bed --dog \
#        --out angsd_bam_trimmed_SE_geno08_transversions_only

################################
#### SHARED, VARIABLE SITES ####
################################

module load EIGENSOFT/7.2.1-gimkl-2018b

in1='angsd_bam_trimmed_SE_geno08_no_damage'

convertf -p <(echo "genotypename:	${in1}.bed
snpname:	        ${in1}.bim
indivname:	      ${in1}.fam
outputformat:	    EIGENSTRAT # important, make sure output format is right 
genotypeoutname:	${in1}.geno
snpoutname:	      ${in1}.snp
indivoutname:	    ${in1}.ind
familynames: NO")

python extract_variable_shared_sites.py

in1='angsd_bam_trimmed_SE_geno08_no_damage_shared_sites'

convertf -p <(echo "genotypename:	${in1}.geno
snpname:	        ${in1}.snp
indivname:	      ${in1}.ind
outputformat:	    PACKEDPED # important, make sure output format is right 
genotypeoutname:	${in1}.bed
snpoutname:	      ${in1}.bim
indivoutname:	    ${in1}.fam")

### SCRIPT: extract_variable_shared_sites.py
#from collections import Counter
#
#def filter_shared_variable_sites(
#    geno_file="angsd_bam_trimmed_SE_geno08_no_damage.geno",
#    snp_file="angsd_bam_trimmed_SE_geno08_no_damage.snp",
#    output_geno="angsd_bam_trimmed_SE_geno08_no_damage_shared_sites.geno",
#    output_snp="angsd_bam_trimmed_SE_geno08_no_damage_shared_sites.snp"
#):
#    variable_indices = []
#    variable_lines = []
#
#    # Step 1: Identify shared variable SNPs
#    with open(geno_file) as f:
#        for idx, line in enumerate(f):
#            line = line.strip()
#            alleles = [g for g in line if g != '9']
#            if len(set(alleles)) > 1:
#                counts = Counter(alleles)
#                sorted_counts = sorted(counts.values(), reverse=True)
#                if len(sorted_counts) >= 2 and sorted_counts[1] >= 2:
#                    variable_indices.append(idx)
#                    variable_lines.append(line)
#
#    # Step 2: Save filtered .geno
#    with open(output_geno, 'w') as g_out:
#        g_out.write('\n'.join(variable_lines) + '\n')
#
#    # Step 3: Save corresponding .snp lines using a set for fast lookup
#    variable_set = set(variable_indices)
#    with open(snp_file) as s_in, open(output_snp, 'w') as s_out:
#        for idx, line in enumerate(s_in):
#            if idx in variable_set:
#                s_out.write(line)
#
#    print(f"Saved {len(variable_lines)} shared variable SNPs.")
#    print(f"  - GENO: {output_geno}")
#    print(f"  - SNP:  {output_snp}")
#
## Run the function
#filter_shared_variable_sites()


########################
### Pre-contact only ###
########################

#Whenua Hou post-contact has the highest missingness. To bump the covered sites, I also created a dataset without these samples.

plink --tfile angsd_bam_trimmed_SE_geno08 --missing-genotype N \
        --make-bed --dog --remove remove_Whenua_Hou_post.txt \
        --out angsd_bam_trimmed_SE_geno08_pre

plink --bfile angsd_bam_trimmed_SE_geno08_pre \
        --make-bed --snps-only --biallelic-only --dog --maf 0.01 \
        --out angsd_bam_trimmed_SE_geno08_pre_filter 
#48,496,657 variants

module load Python/3.11.6-foss-2023a

python identify_CtoT.py
#identified 22,966,033
 
plink --bfile angsd_bam_trimmed_SE_geno08_pre_filter \
        --exclude damage_SE_pre.txt --make-bed --dog \
        --out angsd_bam_trimmed_SE_geno08_pre_no_damage
#25,530,624 variants remaining

#python identify_transitions.py
#
#plink --bfile angsd_bam_trimmed_SE_geno08_pre_filter \
#        --exclude transitions_SE_pre.txt --make-bed --dog \
#        --out angsd_bam_trimmed_SE_geno08_pre_transversions_only