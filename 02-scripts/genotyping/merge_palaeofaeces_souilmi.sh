#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=merge_eigenstrat
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --output merge_eigenstrat.%j.out # CHANGE map1 part each run
#SBATCH --error merge_eigenstrat.%j.err # CHANGE map1 part each runmodule purge 

module load EIGENSOFT/7.2.1-gimkl-2018b

mergeit -p merge_param_file.txt

