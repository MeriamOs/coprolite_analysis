#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=smartPCA
#SBATCH --time=1:00:00
#SBATCH --mem=8GB
#SBATCH --mail-type=ALL
#SBATCH --output smartPCA.%j.out # CHANGE map1 part each run
#SBATCH --error smartPCA.%j.err # CHANGE map1 part each runmodule purge 

module load EIGENSOFT/7.2.1-gimkl-2018b

smartpca -p smartpca.param