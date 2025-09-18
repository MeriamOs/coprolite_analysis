#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J coproid
#SBATCH --time 4:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --output coproid_kuri.%j.out # CHANGE map1 part each run
#SBATCH --error coproid_kuri.%j.err # CHANGE map1 part each runmodule purge

module purge

module load Nextflow/24.10.5

module load Apptainer/1.3.1

export SINGULARITY_TMPDIR=/nesi/nobackup/uoo02328/meriam/container-cache
export SINGULARITY_CACHEDIR=$SINGULARITY_TMPDIR
export NXF_SINGULARITY_CACHEDIR=$SINGULARITY_TMPDIR
setfacl -b "$SINGULARITY_TMPDIR"

#nextflow run nf-core/coproid -r 2.0.0 -profile test,singularity

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/coproid2

nextflow run nf-core/coproid \
	-r 2.0.0 \
	-c coproid.config \
	-profile singularity \
	--outdir '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/coproid2/' \
	--input '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/files/samplesheet2.csv' \
    --genome_sheet "/nesi/nobackup/uoo02328/meriam/coprolite_analysis/files/genomesheet.csv" \
    --sam2lca_db "/home/vanme090/.sam2lca" \
    --kraken2_db "/nesi/nobackup/uoo02328/meriam/kraken/databases/2024-standard" \
	--sp_sources '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/files/sp_sources_september2025.csv' \
	--sp_labels '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/files/sp_labels_september2025.csv' \
	--taxa_sqlite "/home/vanme090/.etetoolkit/taxa.sqlite" \
	--taxa_sqlite_traverse_pkl "/home/vanme090/.etetoolkit/taxa.sqlite.traverse.pkl" \
	--file_prefix 'coproid_kuri' \
    --sam2lca_identity 0.9 \
    --sam2lca_acc2tax 'nucl' \
	-resume
