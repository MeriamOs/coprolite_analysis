#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=smartPCA
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --output smartPCA.%j.out # CHANGE map1 part each run
#SBATCH --error smartPCA.%j.err # CHANGE map1 part each runmodule purge 

module load EIGENSOFT/7.2.1-gimkl-2018b

#smartpca -p smartpca_subset.param
#smartpca -p smartpca_nowolves.param
smartpca -p smartpca_angsd.param