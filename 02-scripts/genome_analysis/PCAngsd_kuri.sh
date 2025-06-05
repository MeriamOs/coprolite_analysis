module purge 
module load Miniconda3

source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

#git clone https://github.com/Rosemeis/pcangsd.git

#cd pcangsd
#conda env create -f environment.yml --solver classic

conda activate pcangsd

#################
### kuri only ###
#################

pcangsd -p pileupcaller.single.maf0.01 -e 2 -t 8 --iter 1000 -o pileupcaller.single.pcangsd 

pcangsd -p angsd_bam_trimmed_SE_geno08_maf0.4_no_damage -e 2 -t 8 --iter 3000 -o angsd_bam_trimmed_SE_geno08_maf0.4_no_damage

pcangsd -p angsd_bam_trimmed_SE_geno08_maf0.33_no_damage -e 2 -t 8 --iter 2000 -o angsd_bam_trimmed_SE_geno08_maf0.33_no_damage

pcangsd -p angsd_bam_trimmed_SE_geno08_maf0.25_no_damage -e 2 -t 8 --iter 3000 -o angsd_bam_trimmed_SE_geno08_maf0.25_no_damage

pcangsd -p angsd_bam_trimmed_SE_geno08_no_damage_shared_sites -e 2 -t 8 --iter 3000 -o angsd_bam_trimmed_SE_geno08_no_damage_shared_sites
