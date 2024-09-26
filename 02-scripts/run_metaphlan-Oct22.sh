#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J metaphlan
#SBATCH --time 3:00:00
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -n 1
#SBATCH --mem=24G
#SBATCH --mail-type=ALL
#SBATCH --output metaphlan_palaeofaeces.%j.out # CHANGE map1 part each run
#SBATCH --error metaphlan_palaeofaeces.%j.err # CHANGE map1 part each runmodule purge

module purge

module load MetaPhlAn/4.0.4-gimkl-2022a-Python-3.10.5

DATADIR='/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/eager/metagenomic_complexity_filter'

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775 ERR5863536 SRR12455959 ERR3761407 ERR3761412;
do

metaphlan ${DATADIR}/${sample}.unmapped.fastq.gz_lowcomplexityremoved.fq.gz \
    --input_type fastq \
    --bowtie2out metaphlan/${sample}_oct22.bt2.out  \
    --nproc 4 \
    --bowtie2db metaphlan/metaphlan_oct22 \
    > metaphlan/${sample}_oct22.metaphlan_profile.txt

done