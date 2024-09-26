#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=eager_Meriam
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --output plink_to_snp.%j.out # CHANGE map1 part each run
#SBATCH --error plink_to_snp.%j.err # CHANGE map1 part each runmodule purge 

module load EIGENSOFT/7.2.1-gimkl-2018b

in1='souilmi2024_bergstrom2022_bergstrom2020_zhang2020_plassais2019_finalset.bergsgrom2022-analyses-loci.tv.c2'

convertf -p <(echo "genotypename:	${in1}.bed
snpname:	        ${in1}.bim
indivname:	      ${in1}.fam
outputformat:	    EIGENSTRAT # important, make sure output format is right 
genotypeoutname:	${in1}.geno
snpoutname:	      ${in1}.snp
indivoutname:	    ${in1}.ind
familynames: NO")

