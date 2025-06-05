module load Python/3.11.6-foss-2023a

python custom_freq_count.py

### SCRIPT: custom_freq_count.py
#import csv
#
#input_vcf = "../angsd_bam_trimmed_SE_no_post_04_haploid.vcf"
#output_frq = "angsd_bam_trimmed_SE_no_post_04_haploid_custom.frq"
#
#with open(input_vcf, "r") as fin, open(output_frq, "w", newline="") as fout:
#    reader = csv.reader(fin, delimiter="\t")
#    writer = csv.writer(fout, delimiter="\t")
#
#    header_written = False
#
#    for row in reader:
#        if row[0].startswith("##"):
#            continue
#        elif row[0] == "#CHROM":
#            sample_names = row[9:]
#            if not header_written:
#                writer.writerow(["CHROM", "POS", "N_ALLELES", "N_CHR", "ALLELE1:FREQ", "ALLELE2:FREQ"])
#                header_written = True
#            continue
#
#        chrom, pos, _, ref, alt = row[0], row[1], row[2], row[3], row[4]
#        genotypes = row[9:]
#
#        ref_count = sum(1 for g in genotypes if g == "0")
#        alt_count = sum(1 for g in genotypes if g == "1")
#        total = ref_count + alt_count
#
#        if total == 0:
#            continue
#
#        ref_freq = ref_count / total
#        alt_freq = alt_count / total
#
#        writer.writerow([
#            chrom,
#            pos,
#            "1",
#            total,
#            f"{ref}:{ref_freq:.3f}",
#            f"{alt}:{alt_freq:.3f}"
#        ])

python calculate_SFS_alt_allele.py 
#Loaded 138,738,634 rows from the .frq file
#Alt Allele Count        Number of Sites
#0	135698128
#1	2913079
#2	117814
#3	6736
#4	2435
#5	404
#6	36
#7	1

### SCRIPT: calculate_SFS_alt_allele.py 
#import pandas as pd
#import matplotlib.pyplot as plt
#
## Load the full file
#df = pd.read_csv("angsd_bam_trimmed_SE_no_post_04_haploid_custom.frq", delim_whitespace=True)
#
#print(f"Loaded {len(df)} rows from the .frq file")
#
## Extract alt allele frequency from ALLELE2:FREQ
#alt_freq = df['ALLELE2:FREQ'].str.split(':').str[1].astype(float)
#n_chr = df['N_CHR'].astype(int)
#
## Calculate alternative allele count
#alt_count = (alt_freq * n_chr).round().astype(int)
#
## Count SFS (Site Frequency Spectrum), including fixed reference (0) and fixed alt (n_chr)
#sfs = alt_count.value_counts().sort_index()
#
#print("\nAlt Allele Count\tNumber of Sites")
#for count, num_snps in sfs.items():
#    print(f"{count}\t{num_snps}")
#
## Plot
#plt.bar(sfs.index, sfs.values, color="purple")
#plt.xlabel("Alternative Allele Count")
#plt.ylabel("Number of Sites")
#plt.title("Site Frequency Spectrum (including fixed sites)")
#plt.tight_layout()
#plt.savefig("angsd_bam_trimmed_SE_no_post_04_haploid_alt_sfs_plot.png", dpi=300)


###########################
#### per arch site SFS ####
###########################

module load BCFtools/1.19-GCC-11.3.0

bcftools view -S ../long_bay_individuals.txt -Ov \
	-o angsd_bam_trimmed_SE_no_post_04_haploid_LB.vcf ../angsd_bam_trimmed_SE_no_post_04_haploid.vcf 

bcftools view -S ../whenua_hou_pre_individuals.txt -Ov \
	-o angsd_bam_trimmed_SE_no_post_04_haploid_WH.vcf ../angsd_bam_trimmed_SE_no_post_04_haploid.vcf 

bcftools view -S ../kahukura_individuals.txt -Ov \
	-o angsd_bam_trimmed_SE_no_post_04_haploid_KH.vcf ../angsd_bam_trimmed_SE_no_post_04_haploid.vcf 

python custom_freq_count.py

python calculate_SFS_alt_allele.py 

###Long Bay
#Loaded 138715240 rows from the .frq file
#
#Alt Allele Count        Number of Sites
#0       137224417
#1       1462242
#2       27788
#3       759
#4       34

### Whenua Hou
#Loaded 138736001 rows from the .frq file
#
#Alt Allele Count        Number of Sites
#0       137567751
#1       1149539
#2       17883
#3       781
#4       46
#5       1

### Kahukuha
#Loaded 121782207 rows from the .frq file
#
#Alt Allele Count        Number of Sites
#0       121312484
#1       466739
#2       2866
#3       113
#4       5