#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=smartPCA
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --output smartPCA.%j.out # CHANGE map1 part each run
#SBATCH --error smartPCA.%j.err # CHANGE map1 part each runmodule purge 

module load EIGENSOFT/7.2.1-gimkl-2018b

nano population_list.txt
#Long_Bay
#Whenua_Hou_post
#Whenua_Hou_pre
#Kahukura

smartpca -p smartpca_maf0.25.param
#genotypename:  angsd_bam_trimmed_SE_geno08_maf0.25_no_damage.bed
#snpname:   angsd_bam_trimmed_SE_geno08_maf0.25_no_damage.bim
#indivname:   angsd_bam_trimmed_SE_geno08_maf0.25_no_damage.fam
#evecoutname:  angsd_bam_trimmed_SE_geno08_maf0.25_no_damage.evec
#evaloutname:   angsd_bam_trimmed_SE_geno08_maf0.25_no_damage.eval
#poplistname: population_list.txt
#lsqproject: NO
#missingmode: NO
#inbreed: YES
#numchrom: 38

smartpca -p smartpca_maf0.33.param
#genotypename:  angsd_bam_trimmed_SE_geno08_maf0.33_no_damage.bed
#snpname:   angsd_bam_trimmed_SE_geno08_maf0.33_no_damage.bim
#indivname:   angsd_bam_trimmed_SE_geno08_maf0.33_no_damage.fam
#evecoutname:  angsd_bam_trimmed_SE_geno08_maf0.33_no_damage.evec
#evaloutname:   angsd_bam_trimmed_SE_geno08_maf0.33_no_damage.eval
#poplistname: population_list.txt
#lsqproject: NO
#missingmode: NO
#inbreed: YES
#numchrom: 38

smartpca -p smartpca_maf0.4.param
#genotypename:  angsd_bam_trimmed_SE_geno08_maf0.4_no_damage.bed
#snpname:   angsd_bam_trimmed_SE_geno08_maf0.4_no_damage.bim
#indivname:   angsd_bam_trimmed_SE_geno08_maf0.4_no_damage.fam
#evecoutname:  angsd_bam_trimmed_SE_geno08_maf0.4_no_damage.evec
#evaloutname:   angsd_bam_trimmed_SE_geno08_maf0.4_no_damage.eval
#poplistname: population_list.txt
#lsqproject: NO
#missingmode: NO
#inbreed: YES
#numchrom: 38

smartpca -p smartpca_shared.param
#genotypename:  angsd_bam_trimmed_SE_geno08_no_damage_shared_sites.bed
#snpname:   angsd_bam_trimmed_SE_geno08_no_damage_shared_sites.bim
#indivname:   angsd_bam_trimmed_SE_geno08_no_damage_shared_sites.fam
#evecoutname:  angsd_bam_trimmed_SE_geno08_no_damage_shared_sites.evec
#evaloutname:   angsd_bam_trimmed_SE_geno08_no_damage_shared_sites.eval
#poplistname: population_list.txt
#lsqproject: NO
#missingmode: NO
#inbreed: YES
#numchrom: 38

smartpca -p smartpca_pileup.param
#genotypename:  pileupcaller.single.maf0.01.bed
#snpname:   pileupcaller.single.maf0.01.bim
#indivname:   pileupcaller.single.maf0.01.fam
#evecoutname:  pileupcaller.single.maf0.01.evec
#evaloutname:   pileupcaller.single.maf0.01.eval
#poplistname: population_list.txt
#lsqproject: NO
#missingmode: NO
#inbreed: YES
#numchrom: 38