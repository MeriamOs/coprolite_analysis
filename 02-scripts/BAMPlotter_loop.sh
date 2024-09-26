ml load Python/3.8.2-gimkl-2020a

ml load SAMtools/1.12-GCC-9.2.0

project_dir=/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/canis_mtDNA

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678
	do
	echo ${sample}
	samtools index ${project_dir}/deduplication/${sample}/${sample}_rmdup.bam  
	python BAMPlotter_modified.py \
    -d ${project_dir}/mapdamage/results_${sample}_rmdup/misincorporation.txt \
	-b ${project_dir}/deduplication/${sample}/${sample}_rmdup.bam \
	-o ${project_dir}/BAMPlotter/${sample}_BAMPlotter_mtDNA.pdf
	done 


#samtools index intermediate_files/${sample}_dedup/${sample}_possort_rmdup.bam  
#-b intermediate_files/${sample}_dedup/${sample}_possort_rmdup.bam
#	-d intermediate_files/MapDamage_results/results_${sample}_possort_rmdup/misincorporation.txt \
#	-b ${sample}_maponly.bam -o ${sample}_BAMPlotter.pdf
