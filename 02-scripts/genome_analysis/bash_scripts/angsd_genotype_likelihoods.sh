#!/bin/bash -e
#SBATCH --cpus-per-task  16
#SBATCH --job-name       glsangsd
#SBATCH --mem            25G
#SBATCH --time           5:00:00
#SBATCH --account        uoo02328
#SBATCH --output         %x_%A_%a.out
#SBATCH --error          %x_%A_%a.err
#SBATCH --hint           nomultithread

module purge
module load angsd/0.935-GCC-9.2.0

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/effective_population

REF="../references/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa"

angsd -bam file_list_angsd.txt -P 32 \
    -GL 2 -doMajorMinor 1 -minMapQ 25 -minQ 20 \
    -doMaf 1 -SNP_pval 1e-4  \
    -doIBS 1 -minind 2 \
    -doCounts 1 -doGlf 2 \
    -uniqueonly 1 -remove_bads 1 \
    -C 50 -baq 1 -makeMatrix 1 -doCov 1 \
    -ref $REF \
    -out angsd_bam_trimmed_GL

module purge
module load SAMtools/1.19-GCC-12.3.0

zcat angsd_bam_trimmed_GL.mafs.gz | cut -f6 |sed 1d >angsd_bam_trimmed_freq

##################
### NGS RELATE ###
##################

./ngsrelate/ngsRelate/ngsRelate -G angsd_bam_trimmed_GL.beagle.gz \
        -n 16 -f angsd_bam_trimmed_freq -O ngsrelate/angsd_bam_trimmed_GL_ngsrelate

module purge

#################
### NGS ADMIX ###
#################

# Define variables
BEAGLE_FILE="angsd_bam_trimmed_GL.beagle.gz"
OUT_PREFIX="angsd_bam_trimmed_GL_ngsAdmix"
THREADS=32
MIN_MAF=0.1                  # Minimum MAF
MIS_TOL=0.8                  # Tolerance for high-quality genotypes
SEEDS=(34545 76490 12321)    # Seeds for reproducibility
MAX_ITER=2000                # Maximum number of EM iterations

# Loop through K values
for K in {2..8}; do
    echo "Running ngsAdmix for K=${K}"
    
    # Loop through seeds
    for SEED in "${SEEDS[@]}"; do
        echo "Using seed=${SEED}"
        
        # Run ngsAdmix
        ./ngsadmix/NGSadmix/NGSadmix -likes $BEAGLE_FILE \
                 -K $K \
                 -outfiles ngsadmix/"${OUT_PREFIX}_K${K}_seed${SEED}" \
                 -seed $SEED \
                 -minMaf $MIN_MAF \
                 -misTol $MIS_TOL \
                 -P $THREADS \
                 -maxiter $MAX_ITER
    done
done

#################
### PCA ANGSD ###
#################

module purge 
module load Miniconda3

source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

#git clone https://github.com/Rosemeis/pcangsd.git

#cd pcangsd
#conda env create -f environment.yml --solver classic

conda activate pcangsd

### Kuri genotype likelihoods ###

pcangsd -b angsd_bam_trimmed_GL.beagle.gz --inbreed-sites --inbreed-samples \
    -e 2 -t 8 --iter 1000 -o angsd_bam_trimmed_GL
