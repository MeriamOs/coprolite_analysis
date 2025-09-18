#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=contammix_Meriam
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --output contammix.%j.out # CHANGE map1 part each run
#SBATCH --error contammix.%j.err # CHANGE map1 part each runmodule purge #prepare contaminants file

#create multiple fasta file contaminants.fasta with pig reference (EF545567), human (Homo sapiens, NC_012920.1), 
#cow (Bos taurus, Gen Bank NC_006853.1), dog (Canis familiaris, NC_002008.4), and chicken (Gallus gallus, Gen Bank NC_001323.1)

module purge
module load MAFFT/7.505-gimkl-2022a-with-extensions
module load R/4.2.1-gimkl-2022a
module load SAMtools/1.16.1-GCC-11.3.0
module load BWA/0.7.17-GCC-11.3.0

mafft contaminants.fasta > contaminants_aligned.fasta

BAMDIR=/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/canis_mtDNA/deduplication
CONDIR=/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/canis_mtDNA/consensus_sequence

for SAMPLE in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
            MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
do 
    cat ${CONDIR}/${SAMPLE}.fasta contaminants_aligned.fasta > ${SAMPLE}_contaminants.fasta
    mafft ${SAMPLE}_contaminants.fasta > ${SAMPLE}_contaminants_aligned.fasta
    bwa index ${CONDIR}/${SAMPLE}.fasta 
    samtools fastq ${BAMDIR}/${SAMPLE}/${SAMPLE}_rmdup.bam > ${SAMPLE}_rmdup.fq
    bwa aln -t 24 -n 0.03 -l 1024 ${CONDIR}/${SAMPLE}.fasta ${SAMPLE}_rmdup.fq > ${SAMPLE}.sai
    bwa samse ${CONDIR}/${SAMPLE}.fasta ${SAMPLE}.sai ${SAMPLE}_rmdup.fq -f ${SAMPLE}.sam
    samtools view -b -S -q 20 ${SAMPLE}.sam > ${SAMPLE}_contammix.bam
    contamMix/exec/estimate.R --samFn ${SAMPLE}_contammix.bam \
    --malnFn ${SAMPLE}_contaminants_aligned.fasta --figure ${SAMPLE}_contamMix \
    --nrThreads 24 \
    > "${SAMPLE}_log.txt" 2>&1
done 

for file in *_log.txt; 
    do 
    name="${file%_log.txt}"
    echo -n "$name: "
    grep 'MAP' $file | awk '{print $3}'
done
