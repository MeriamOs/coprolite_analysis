#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=eager_canFam3.1
#SBATCH --time=168:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --output eager_canFam3.1.%j.out # CHANGE map1 part each run
#SBATCH --error eager_canFam3.1.%j.err # CHANGE map1 part each runmodule purge 

module purge

module load Nextflow/22.04.3

module load Singularity/3.11.3

export SINGULARITY_TMPDIR=/nesi/nobackup/uoo02328/meriam/container-cache
export SINGULARITY_CACHEDIR=$SINGULARITY_TMPDIR
export NXF_SINGULARITY_CACHEDIR=$SINGULARITY_TMPDIR
setfacl -b "$SINGULARITY_TMPDIR"

#nextflow run nf-core/eager -profile test,singularity

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis

#nextflow run nf-core/eager \
#-r 2.5.0 \
#-c 02_eager.config \
#-profile singularity \
#--outdir '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/eager/' \
#--input '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/03-data/whenua_hou/Blank*_R{1,2}_001.fastq.gz' \
#--fasta '/nesi/project/uoo02328/references/domesticates_mtDNA/Canis_familiaris_genome_GCF_014441545.1.fna' \
#--skip_adapterremoval --skip_preseq \
#--bwaalnn 0.03 --dedupper dedup --mergedonly \
#--damage_calculation_tool 'mapdamage' --mapdamage_downsample 100000 \
#--run_bam_filtering --bam_unmapped_type 'fastq' --metagenomic_complexity_filter \
#-resume 

nextflow run nf-core/eager \
-r 2.5.0 \
-c 02_eager.config \
-profile singularity \
--outdir '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/eager_CanFam3.1/' \
--input '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/eager/adapterremoval/output/MS*.fq.gz' \
--fasta '/nesi/project/uoo02328/references/domesticates_mtDNA/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa' \
--bwa_index '/nesi/project/uoo02328/references/domesticates_mtDNA/' \
--save_reference \
--single_end \
--skip_adapterremoval --skip_preseq \
--skip_damage_calculation \
--run_trim_bam --bamutils_clip_single_stranded_none_udg_left 2 \
--bamutils_clip_single_stranded_none_udg_left 2 \
--bamutils_softclip \
--bwaalnn 0.03 --dedupper dedup --mergedonly \
--run_genotyping --genotyping_tool 'pileupcaller' \
--genotyping_source 'trimmed' \
--pileupcaller_bedfile '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/genotyping/souilmi2024_bergstrom2022_bergstrom2020_zhang2020_plassais2019_finalset.bergsgrom2022-analyses-loci.convert.bed' \
--pileupcaller_snpfile '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/genotyping/souilmi2024_bergstrom2022_bergstrom2020_zhang2020_plassais2019_finalset.bergsgrom2022-analyses-loci.converted.snp' \
--pileupcaller_method 'randomHaploid' \
--run_bcftools_stats \
-resume 

#--pileupcaller_bedfile '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/genotyping/canFam3.1.bed' \

#--run_genotyping --genotyping_tool 'ug' \
#--gatk_ploidy 1 --gatk_ug_out_mode 'EMIT_ALL_SITES' --gatk_ug_genotype_model 'SNP' \
#--run_vcf2genome --vcf2genome_minc 3 --vcf2genome_minfreq 0.66 \

#clean up intermediate files
nextflow clean -f -k

