cp ../souilmi_kuri_subset_to_keep.txt dingo_ancient.txt
#Initially I kept the ancient dingoes in list - but ended up romoving these because got an error
#Error in nj(mdist.dist) : 
#  missing values are not allowed in the distance matrix
#Consider using njs()
#Execution halted

NAME="merged_kuri_souilmi"

plink --bfile ../${NAME} --make-bed --remove dingo_ancient.txt --dog \
	--out merged_kuri_souilmi_ibs

#module load EIGENSOFT/7.2.1-gimkl-2018b
module purge
module load PLINK/1.09b6.16

#plink --distance square 1-ibs --out merged_kuri_souilmi_ibs_matrix --bfile ../merged_kuri_souilmi --dog

#plink --distance square 1-ibs --bfile angsd_bam_trimmed_geno02_transversions_only --out angsd_bam_transversions_ibs_matrix --dog

#plink --distance square 1-ibs --bfile ../genotyping/pileupcaller_analysis/variable_sites --out ibs_matrix/pileupcaller_variable_ibs_matrix --dog

#############################
### GLOBAL DOGS BOOTSTRAP ###
#############################

NAME="merged_kuri_souilmi_ibs"
OUT="bootstrap"

module load R/4.3.2-foss-2023a

plink --bfile ${NAME} --distance square 1-ibs --dog --out ${OUT}/${NAME}
cut -f 2 ${NAME}.bim > ${OUT}/site_list.txt
cut -f 2 ${OUT}/${NAME}.mdist.id |  perl -pe 's/\n/\t/g' | awk 'BEGIN{FS=OFS="ID\t"}{print value OFS $0}' > ${OUT}/matrix_head #Converts list to a tab-delimited line, and adds an 'ID' column at start (to account for offset)
cut -f 2 ${OUT}/${NAME}.mdist.id > ${OUT}/matrix_labels
paste ${OUT}/matrix_labels ${OUT}/${NAME}.mdist | cat ${OUT}/matrix_head - > ${OUT}/${NAME}.d

Rscript bionj.R ${OUT}/${NAME}.d ${OUT}/${NAME}.tree

#bionj.R
	#args <- commandArgs(trailingOnly=TRUE)
	#library(ape)
	#distm<-read.table(args[1], header = TRUE, row.names = 1)
	#mdist.dist <- as.dist(distm)
	#nj.tree <- nj(mdist.dist)
	#max.dist <- max(node.depth.edgelength(nj.tree))
	#root.node <- which(node.depth.edgelength(nj.tree) == max.dist)
	#rooted.tree <- root(nj.tree, root.node)
	#write.tree(rooted.tree,file=args[2])

for i in {1..1000}; 
do
    #Select 50,000 or 100,000 random sites to use to generate each tree
	echo -e "\n Subsetting Sites for Tree ${i}";
	cat ${OUT}/site_list.txt | shuf | head -100000 > trees/sites_rep_${i}.txt; 
	
	#Calculate IBS distances between each sample in reduced dataset, and removes sites file afterwards
	echo -e "\n Creating IBS Distance Matrix for Tree ${i}";
	plink --bfile ${NAME} --distance square 1-ibs --extract trees/sites_rep_${i}.txt --dog --allow-no-sex --make-bed --out trees/${NAME}_rep_${i};
	rm trees/sites_rep_${i}.txt;
	
	#Formatting PLINK output
	cut -f 2 trees/${NAME}_rep_${i}.mdist.id |  perl -pe 's/\n/\t/g' | awk 'BEGIN{FS=OFS="ID\t"}{print value OFS $0}' > trees/matrix_head; #Converts list to a tab-delimted line, and adds an 'ID' column at start (to account for offset)
	cut -f 2 trees/${NAME}_rep_${i}.mdist.id > trees/matrix_labels;
	paste trees/matrix_labels trees/${NAME}_rep_${i}.mdist | cat trees/matrix_head - > trees/${NAME}_rep_${i}.d;

	#Create a tree specifying input/output files
	echo -e "\n Creating Tree ${i}";
	Rscript bionj.R trees/${NAME}_rep_${i}.d trees/${NAME}_rep_${i}.tree;

	#Create consensus tree with bootstrap
	cat trees/${NAME}_rep_${i}.tree | head -1 >> ${OUT}/${NAME}_bootstrap.tree;
	echo -e "\n Finished Tree ${i}";
done

module purge 
module load Miniconda3

source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

#conda create -n booster booster=0.1.2 --solver=classic
conda activate booster

cat trees/${NAME}_rep*.tree > ${OUT}/all_replicate.tree
booster -i ${OUT}/${NAME}_bootstrap.tree -b ${OUT}/all_replicate.tree -o ${OUT}/${NAME}_BOOSTER_tree.nwk


######################
### KURI BOOTSTRAP ###
######################

plink --bfile ../pileupcaller.single --make-bed --dog --out pileupcaller.single.maf0.01 --maf 0.01 \
	--snps-only --biallelic-only --missing-genotype N --geno 0.8
#60,864 variants remain

plink --bfile ../../effective_population/angsd_bam_trimmed_SE_geno08_no_damage --make-bed --dog \
	--out angsd_bam_trimmed_SE_geno04_no_damage \
	--snps-only --biallelic-only --missing-genotype N --geno 0.4
#529,549 variants remain
#there were only 5185 variants for --geno 0.2
#24,023 for --geno 0.3

#there was no structure when using the prefiltered --maf 0.01 setting, so filtered further to only
#retain shared SNPs
plink --bfile ../../effective_population/angsd_bam_trimmed_SE_geno08_no_damage --make-bed --dog \
	--out angsd_bam_trimmed_SE_geno08_maf0.25_no_damage --maf 0.25 --snps-only --biallelic-only \
	--missing-genotype N --geno 0.8
#5,135,835 remain

plink --bfile ../../effective_population/angsd_bam_trimmed_SE_geno08_no_damage --make-bed --dog \
	--out angsd_bam_trimmed_SE_geno08_maf0.33_no_damage --maf 0.33 --snps-only --biallelic-only \
	--missing-genotype N --geno 0.8
#486,631 remain

plink --bfile ../../effective_population/angsd_bam_trimmed_SE_geno08_no_damage --make-bed --dog \
	--out angsd_bam_trimmed_SE_geno08_maf0.4_no_damage --maf 0.4 --snps-only --biallelic-only \
	--missing-genotype N --geno 0.8
#257,900

#Run the below twice, once for the pieupcaller dataset, and the variously filtered angsd no damage dataset
#NAME="pileupcaller.single.maf0.01"
#NAME="angsd_bam_trimmed_SE_geno08_maf0.25_no_damage"
#NAME="angsd_bam_trimmed_SE_geno08_maf0.33_no_damage"
#NAME="angsd_bam_trimmed_SE_geno08_maf0.4_no_damage"
NAME="angsd_bam_trimmed_SE_geno08_no_damage_shared_sites"
OUT="kuri_only"

module load R/4.3.2-foss-2023a

plink --bfile ${NAME} --distance square 1-ibs --dog --out ${OUT}/${NAME}
cut -f 2 ${NAME}.bim > ${OUT}/site_list.txt
cut -f 2 ${OUT}/${NAME}.mdist.id |  perl -pe 's/\n/\t/g' | awk 'BEGIN{FS=OFS="ID\t"}{print value OFS $0}' > ${OUT}/matrix_head #Converts list to a tab-delimted line, and adds an 'ID' column at start (to account for offset)
cut -f 2 ${OUT}/${NAME}.mdist.id > ${OUT}/matrix_labels
paste ${OUT}/matrix_labels ${OUT}/${NAME}.mdist | cat ${OUT}/matrix_head - > ${OUT}/${NAME}.d

Rscript ${OUT}/bionj.R ${OUT}/${NAME}.d ${OUT}/${NAME}.tree

#bionj.R
    #args <- commandArgs(trailingOnly=TRUE)
    #library(ape)
    #library(phangorn)
    #distm<-read.table(args[1], header = TRUE, row.names = 1)
    #mdist.dist <- as.dist(distm)
    #nj.tree <- nj(mdist.dist)
    #rooted.tree <- midpoint(nj.tree)
    #write.tree(rooted.tree,file=args[2])

for i in {1..1000}; 
do
    #Select 50,000 or 100,000 random sites to use to generate each tree
	echo -e "\n Subsetting Sites for Tree ${i}";
	cat ${OUT}/site_list.txt | shuf | head -100000 > ${OUT}/trees/sites_rep_${i}.txt; 
	
	#Calculate IBS distances between each sample in reduced dataset, and removes sites file afterwards
	echo -e "\n Creating IBS Distance Matrix for Tree ${i}";
	plink --bfile ${NAME} --distance square 1-ibs --extract ${OUT}/trees/sites_rep_${i}.txt --dog --allow-no-sex --make-bed --out ${OUT}/trees/${NAME}_rep_${i};
	rm ${OUT}/trees/sites_rep_${i}.txt;
	
	#Formatting PLINK output
	cut -f 2 ${OUT}/trees/${NAME}_rep_${i}.mdist.id |  perl -pe 's/\n/\t/g' | awk 'BEGIN{FS=OFS="ID\t"}{print value OFS $0}' > ${OUT}/trees/matrix_head; #Converts list to a tab-delimted line, and adds an 'ID' column at start (to account for offset)
	cut -f 2 ${OUT}/trees/${NAME}_rep_${i}.mdist.id > ${OUT}/trees/matrix_labels;
	paste ${OUT}/trees/matrix_labels ${OUT}/trees/${NAME}_rep_${i}.mdist | cat ${OUT}/trees/matrix_head - > ${OUT}/trees/${NAME}_rep_${i}.d;

	#Create a tree specifying input/output files
	echo -e "\n Creating Tree ${i}";
	Rscript ${OUT}/bionj.R ${OUT}/trees/${NAME}_rep_${i}.d ${OUT}/trees/${NAME}_rep_${i}.tree;

	#Create consensus tree with bootstrap
	cat ${OUT}/trees/${NAME}_rep_${i}.tree | head -1 >> ${OUT}/${NAME}_bootstrap.tree;
	echo -e "\n Finished Tree ${i}";
done

module purge 
module load Miniconda3

source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

#conda create -n booster booster=0.1.2 --solver=classic
conda activate booster

#Add bootstrap values to tree
cat ${OUT}/trees/${NAME}_rep*.tree > ${OUT}/all_replicate.tree
booster -i ${OUT}/${NAME}_bootstrap.tree -b ${OUT}/all_replicate.tree \
	-o ${OUT}/${NAME}_BOOSTER_tree.nwk