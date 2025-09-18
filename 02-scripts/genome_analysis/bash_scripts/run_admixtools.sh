module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

#conda create -n admixtools -c bioconda -c conda-forge admixtools --solver classic

conda activate admixtools

ln -s ../genotyping/merged_kuri_souilmi.* .

module load EIGENSOFT/7.2.1-gimkl-2018b

in1='merged_kuri_souilmi'

convertf -p <(echo "
genotypename:	    ${in1}.bed
snpname:	        ${in1}.bim
indivname:	        ${in1}.fam
outputformat:	    EIGENSTRAT # important, make sure output format is right 
genotypeoutname:    ${in1}.geno
snpoutname:	        ${in1}.snp
indivoutname:	    ${in1}.ind
familynames:        NO")

#accidentily ran with familynames: NO, so manually added the family names to the ind file
awk 'NR==FNR {a[NR]=$1; next} { $3 = a[FNR]; print }' merged_kuri_souilmi.fam merged_kuri_souilmi.ind \
    > merged_kuri_souilmi_modified.ind

qp3Pop -p qp3pop.params >shared_drift_qp3pop.logfile

### FILE: qp3pop.params
#genotypename:  merged_kuri_souilmi.geno
#snpname:   merged_kuri_souilmi.snp
#indivname:   merged_kuri_souilmi_modified.ind
#popfilename:  list_qp3test.txt
#inbreed: YES
#numchrom: 38
#outgroupmode: NO
#allsnps: YES

qpDstat -p qpDstat_European_breed.params >admix_qpDstat_european.logfile
#D(Jackal, Post-contact Kuri; Pre-contact Kuri, European Breed)

### FILE: qpDstat_European_breed.params
#genotypename:  merged_kuri_souilmi.geno
#snpname:   merged_kuri_souilmi.snp
#indivname:   merged_kuri_souilmi_samples.ind
#popfilename:  list_qpDstat_European.txt
#inbreed: YES
#numchrom: 38
#outgroupmode: YES

qpDstat -p qpDstat_wolves.params >admix_qpDstat_wolves.logfile
#D(Jackal, Wolf population, Kuri, European Breed)

### FILE: qpDstat_wolves.params
#genotypename:  merged_kuri_souilmi.geno
#snpname:   merged_kuri_souilmi.snp
#indivname:   merged_kuri_souilmi_modified.ind
#popfilename:  list_qpDstat_wolves.txt
#inbreed: YES
#numchrom: 38
#outgroupmode: YES

