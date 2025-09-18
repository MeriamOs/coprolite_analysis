#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=READ
#SBATCH --time=04:00:00
#SBATCH --mem=8GB
#SBATCH --output READ.%j.out # CHANGE map1 part each run
#SBATCH --error READ.%j.err # CHANGE map1 part each runmodule purge

module purge 
module load Miniconda3

source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

#conda create -n readv2 python=3 pandas numpy pip matplotlib 

conda activate readv2

#conda install plinkio

#module load PLINK/2.00a2.3
plink2 --bfile pileupcaller.double --maf 0.01 --dog --make-bed --out pileupcaller.double_maf0.01
#--maf minimum allele frequency

python READv2/READ2.py -i pileupcaller.double_maf0.01 
mv meansP0_AncientDNA_normalized_READv2 meansP0_AncientDNA_normalized_maf0.01_median
mv Read_Results.tsv Read_Results_maf0.01_median.tsv
mv READ.pdf READ_maf0.01_median.pdf

#-n/--norm_method is the normalisation method: when using median (default) it assumes most pairs are unrelated. 
    #Mean uses the mean across all pairs, which is more sensitive to outliers (e.g. identical twins) - tested as some might be from same individual
python READv2/READ2.py -i pileupcaller.double_maf0.01 -n mean
mv meansP0_AncientDNA_normalized_READv2 meansP0_AncientDNA_normalized_maf0.01_mean
mv Read_Results.tsv Read_Results_maf0.01_mean.tsv
mv READ.pdf READ_maf0.01_mean.pdf

##focus on results with at least 10,000 overlapping SNPs