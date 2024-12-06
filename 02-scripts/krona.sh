module purge

module load MetaPhlAn/4.1.0-gimkl-2022a-Python-3.10.5

for sample in MS11684;
do

python /opt/nesi/CS400_centos7_bdw/MetaPhlAn/4.1.0-gimkl-2022a-Python-3.10.5/bin/metaphlan2krona.py -p ${sample}_4.1.0.metaphlan_profile.txt -k krona/${sample}_krona.out

ktImportText -o krona/${sample}_krona.html krona/${sample}_krona.out

done