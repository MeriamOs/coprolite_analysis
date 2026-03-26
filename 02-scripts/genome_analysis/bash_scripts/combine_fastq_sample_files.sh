## NEXTSEQ DATA

mkdir nextseq
cd nextseq

for ID in 680636061 680636062 680636063 680636064 680636065 680636066 680636067 680636068 \
        686504976 686504977 686504978 686504979 686504980 686504981 686504982 686504983 \
        680636069 680636070 686504984 686504985 686504986 686504987 #blanks
    do 
    /nesi/project/uoo02328/programs/bs download biosample  -i ${ID} --extension=fastq.gz
done

for file in *.gz; do mv "$file" "${file/-00-01_S*_L001/}"; done
for file in *.gz; do mv "$file" "${file/AACJ27FHV-/}"; done
for file in *.gz; do mv "$file" "${file/AACJ7KGHV-/}"; done
for file in *.gz; do mv "$file" "${file/_001}"; done

cd ..

## NANORUN DATA 

mkdir nanoruns
cd nanoruns

for ID in 644371793 644371794 644371797 644371798 644371799 644371800 644371801 644371802 \
        644371803 644371807 644371808 644371810 676215598 676215599 676215600 676215601 \
        676215602 676215603 676215604 676215605 680636635 680636636 680636639 680636640 \
        684927845 684927846 684927847 684927848 684927849 684927850 684927851 684927852 \
        676215606 676215607 684927853 684927854 684927855 684927856 #blanks
    do 
    /nesi/project/uoo02328/programs/bs download biosample  -i ${ID} --extension=fastq.gz
done

for file in *.gz; do echo mv "$file" "${file/-00-01_S*_L001/}"; done

cd ..

cd ..

mkdir renamed
mv nanoruns/*.gz renamed
mv nextseq/*.gz renamed
cd renamed

## Whenua Hou

cat	8511-51_R1_001.fastq.gz	8889-23_R1_001.fastq.gz	8889-23_R1.fastq.gz	>	MS11669_R1_001.fastq.gz
cat	8511-51_R2_001.fastq.gz	8889-23_R2_001.fastq.gz	8889-23_R2.fastq.gz	>	MS11669_R2_001.fastq.gz
cat	8511-52_R1_001.fastq.gz	8889-24_R1_001.fastq.gz	8889-24_R1.fastq.gz	>	MS11670_R1_001.fastq.gz
cat	8511-52_R2_001.fastq.gz	8889-24_R2_001.fastq.gz	8889-24_R2.fastq.gz	>	MS11670_R2_001.fastq.gz
cat	8511-55_R1_001.fastq.gz	8889-25_R1_001.fastq.gz	8889-25_R1.fastq.gz	>	MS11673_R1_001.fastq.gz
cat	8511-55_R2_001.fastq.gz	8889-25_R2_001.fastq.gz	8889-25_R2.fastq.gz	>	MS11673_R2_001.fastq.gz
cat	8511-56_R1_001.fastq.gz	8889-26_R1_001.fastq.gz	8889-26_R1.fastq.gz	>	MS11674_R1_001.fastq.gz
cat	8511-56_R2_001.fastq.gz	8889-26_R2_001.fastq.gz	8889-26_R2.fastq.gz	>	MS11674_R2_001.fastq.gz
cat	8511-57_R1_001.fastq.gz	8889-27_R1_001.fastq.gz	8889-27_R1.fastq.gz	>	MS11675_R1_001.fastq.gz
cat	8511-57_R2_001.fastq.gz	8889-27_R2_001.fastq.gz	8889-27_R2.fastq.gz	>	MS11675_R2_001.fastq.gz
cat	8511-58_R1_001.fastq.gz	8889-28_R1_001.fastq.gz	8889-28_R1.fastq.gz	>	MS11676_R1_001.fastq.gz
cat	8511-58_R2_001.fastq.gz	8889-28_R2_001.fastq.gz	8889-28_R2.fastq.gz	>	MS11676_R2_001.fastq.gz
cat	8511-59_R1_001.fastq.gz	8889-29_R1_001.fastq.gz	8889-29_R1.fastq.gz	>	MS11677_R1_001.fastq.gz
cat	8511-59_R2_001.fastq.gz	8889-29_R2_001.fastq.gz	8889-29_R2.fastq.gz	>	MS11677_R2_001.fastq.gz
cat	8511-60_R1_001.fastq.gz	8889-30_R1_001.fastq.gz	8889-30_R1.fastq.gz	>	MS11678_R1_001.fastq.gz
cat	8511-60_R2_001.fastq.gz	8889-30_R2_001.fastq.gz	8889-30_R2.fastq.gz	>	MS11678_R2_001.fastq.gz

## Long Bay

cat	8511-61_R1_001.fastq.gz	8842-21_R1_001.fastq.gz	8842-21_R1.fastq.gz	>	MS11679_R1_001.fastq.gz
cat	8511-61_R2_001.fastq.gz	8842-21_R2_001.fastq.gz	8842-21_R2.fastq.gz	>	MS11679_R2_001.fastq.gz
cat	8511-65_R1_001.fastq.gz	8842-22_R1_001.fastq.gz	8842-22_R1.fastq.gz	>	MS11683_R1_001.fastq.gz
cat	8511-65_R2_001.fastq.gz	8842-22_R2_001.fastq.gz	8842-22_R2.fastq.gz	>	MS11683_R2_001.fastq.gz
cat	8511-66_R1_001.fastq.gz	8842-23_R1_001.fastq.gz	8842-23_R1.fastq.gz	>	MS11684_R1_001.fastq.gz
cat	8511-66_R2_001.fastq.gz	8842-23_R2_001.fastq.gz	8842-23_R2.fastq.gz	>	MS11684_R2_001.fastq.gz
cat	8511-68_R1_001.fastq.gz	8842-24_R1_001.fastq.gz	8842-24_R1.fastq.gz	>	MS11686_R1_001.fastq.gz
cat	8511-68_R2_001.fastq.gz	8842-24_R2_001.fastq.gz	8842-24_R2.fastq.gz	>	MS11686_R2_001.fastq.gz

## Kahukura 

cat	8842-12_R1_001.fastq.gz	8842-25_R1_001.fastq.gz	8842-25_R1.fastq.gz	>	MS11770_R1_001.fastq.gz
cat	8842-12_R2_001.fastq.gz	8842-25_R2_001.fastq.gz	8842-25_R2.fastq.gz	>	MS11770_R2_001.fastq.gz
cat	8842-13_R1_001.fastq.gz	8842-26_R1_001.fastq.gz	8842-26_R1.fastq.gz	>	MS11771_R1_001.fastq.gz
cat	8842-13_R2_001.fastq.gz	8842-26_R2_001.fastq.gz	8842-26_R2.fastq.gz	>	MS11771_R2_001.fastq.gz
cat	8842-16_R1_001.fastq.gz	8842-27_R1_001.fastq.gz	8842-27_R1.fastq.gz	>	MS11774_R1_001.fastq.gz
cat	8842-16_R2_001.fastq.gz	8842-27_R2_001.fastq.gz	8842-27_R2.fastq.gz	>	MS11774_R2_001.fastq.gz
cat	8842-17_R1_001.fastq.gz	8842-28_R1_001.fastq.gz	8842-28_R1.fastq.gz	>	MS11775_R1_001.fastq.gz
cat	8842-17_R2_001.fastq.gz	8842-28_R2_001.fastq.gz	8842-28_R2.fastq.gz	>	MS11775_R2_001.fastq.gz

## BLANKS

cat	8889-31_R1_001.fastq.gz	8889-31_R1.fastq.gz	>	Blank1_WH_R1_001.fastq.gz
cat	8889-31_R2_001.fastq.gz	8889-31_R2.fastq.gz	>	Blank1_WH_R2_001.fastq.gz
cat	8889-32_R1_001.fastq.gz	8889-32_R1.fastq.gz	>	Blank2_WH_R1_001.fastq.gz
cat	8889-32_R2_001.fastq.gz	8889-32_R2.fastq.gz	>	Blank2_WH_R2_001.fastq.gz
cat	8842-29_R1_001.fastq.gz	8842-29_R1.fastq.gz	>	LB_blank_1_R1_001.fastq.gz
cat	8842-29_R2_001.fastq.gz	8842-29_R2.fastq.gz	>	LB_blank_1_R2_001.fastq.gz
cat	8842-30_R1_001.fastq.gz	8842-30_R1.fastq.gz	>	LB_blank_2_R1_001.fastq.gz
cat	8842-30_R2_001.fastq.gz	8842-30_R2.fastq.gz	>	LB_blank_2_R2_001.fastq.gz
cat	8842-31_R1_001.fastq.gz	8842-31_R1.fastq.gz	>	KH_blank_1_R1_001.fastq.gz
cat	8842-31_R2_001.fastq.gz	8842-31_R2.fastq.gz	>	KH_blank_1_R2_001.fastq.gz
cat	8842-32_R1_001.fastq.gz	8842-32_R1.fastq.gz	>	KH_blank_2_R1_001.fastq.gz
cat	8842-32_R2_001.fastq.gz	8842-32_R2.fastq.gz	>	KH_blank_2_R2_001.fastq.gz
