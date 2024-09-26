#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J krk2
#SBATCH --time 8:00:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -n 1
#SBATCH --mem=200G
#SBATCH --mail-type=ALL
#SBATCH --output krk2_palaeofaeces.%j.out # CHANGE map1 part each run
#SBATCH --error krk2_palaeofaeces.%j.err # CHANGE map1 part each run

module purge

module load Kraken2
# note: adding blast module to get dustmasker program to finish
module load BLAST/2.13.0-GCC-11.3.0
## go to directory
## make the following folder first: 
## mkdir /nesi/nobackup/uoo02328/meriam/kraken/all_reads

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis

DATADIR='/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/eager/metagenomic_complexity_filter'
#'/nesi/nobackup/uoo02328/meriam/coprolites/comparative_data'
#'/nesi/project/uoo02328/meriam/coprolite_data/results'
#'/nesi/project/uoo02328/meriam/pathogen_data/HiSeq_data_OG7401/results'
#'/nesi/nobackup/uoo02328/meriam/coprolites'

## reference
ref='/nesi/nobackup/uoo02328/meriam/kraken/databases/2024-palaeofaeces'
ref2='/opt/nesi/db/Kraken2/nt'
#'/nesi/nobackup/uoo02328/meriam/kraken/databases/2024-palaeofaeces'
#'/nesi/nobackup/uoo02328/meriam/kraken/databases/2023-plusPFP-16GB'
#'/opt/nesi/db/Kraken2/nt'
#'/nesi/nobackup/uoo02328/meriam/kraken/databases/2022-palaeofaeces'
#'/nesi/nobackup/uoo02328/meriam/kraken/databases/2021-plusPFP'
#'/nesi/nobackup/uoo02328/meriam/kraken/databases/2020-standard'
#'/nesi/nobackup/uoo02328/meriam/kraken'
#'/opt/nesi/db/Kraken2/standard-2018-09'
#ref='/nesi/project/uoo02328/programs/minikraken2_v2_8GB_database/minikraken2_v2_8GB_201904'
#ref='/PATH/TO/8GB_DATABASE/'

## run kraken

#MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 Blank1_WH Blank2_WH
#MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775 KH_blank_1 KH_blank_2 LB_blank_1 LB_blank_2
#ERR5863536 SRR12455959 ERR3761407 ERR3761412

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775 ERR5863536 SRR12455959 ERR3761407 ERR3761412;
do
kraken2 --db $ref2 \
	--threads 32 \
	--classified-out kraken-filtered/${sample}_kraken_contig_conf0.50.fasta \
	--report kraken-filtered/${sample}_kraken_conf0.50.txt \
	--output kraken-filtered/${sample}_kraken_conf0.50 \
	--use-names \
	--gzip-compressed \
	--minimum-base-quality 30 \
	--confidence 0.50 \
	${DATADIR}/${sample}.unmapped.fastq.gz_lowcomplexityremoved.fq.gz 
done
    
#for sample in MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775 KH_blank_1 KH_blank_2 LB_blank_1 LB_blank_2;
#do

#kraken2 --db $ref2 \
#	--threads 32 \
#	--classified-out kraken-nt/${sample}_kraken_contig_conf0.50.fasta \
#	--report kraken-nt/${sample}_kraken_conf0.50.txt \
#	--output kraken-nt/${sample}_kraken_conf0.50 \
#	--use-names \
#	--gzip-compressed \
#	--minimum-base-quality 30 \
#	--confidence 0.50 \
#	${DATADIR}/${sample}_R1_001.fastq_L0.pe.combined.fq.gz
#done