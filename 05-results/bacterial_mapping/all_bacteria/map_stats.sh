ref1='AP019711.1'
refname1='Amedibacterium_intestinale'
ref2='CP075979.1'
refname2='Clostridium_perfringens'
ref3='U00096.3'
refname3='Escherichia_coli'
ref4='AP012344.1'
refname4='Helicobacter_cinaedi'
ref5='AP024814.1'
refname5='Helicobacter_NHP19-003'
ref6='AP025240.1'
refname6='Phocaeicola_vulgatus'
ref7='CP023536.1'
refname7='Providencia_alcalifaciens'
ref8='AP019711.1'
refname8='Amedibacterium_intestinale'
ref9='NZ_CP118962.1'
refname9='Enterococcus_faecalis'
ref10='NZ_CP038996.1'
refname10='Enterococcus_faecium'
ref11='NZ_CP013119.1'
refname11='Alcaligenes_faecalis'
ref12='NZ_AP031449.1'
refname12='Bacteroides_stercoris'
ref13='NZ_CP017279.1'
refname13='Enterobacter_ludwigii'
ref14='NZ_CP101467.1'
refname14='Rhodococcus_gordoniae'
ref15='CP094882.1'
refname15='Lactococcus_petauri'

samtools coverage MS11669_maponly.bam | head -n 1 > header.txt

cat header.txt | sed "s/#rname/$refname1/" > ${refname1}_stats.txt
cat header.txt | sed "s/#rname/$refname2/" > ${refname2}_stats.txt
cat header.txt | sed "s/#rname/$refname3/" > ${refname3}_stats.txt
cat header.txt | sed "s/#rname/$refname4/" > ${refname4}_stats.txt
cat header.txt | sed "s/#rname/$refname5/" > ${refname5}_stats.txt
cat header.txt | sed "s/#rname/$refname6/" > ${refname6}_stats.txt
cat header.txt | sed "s/#rname/$refname7/" > ${refname7}_stats.txt
cat header.txt | sed "s/#rname/$refname8/" > ${refname8}_stats.txt
cat header.txt | sed "s/#rname/$refname9/" > ${refname9}_stats.txt
cat header.txt | sed "s/#rname/$refname10/" > ${refname10}_stats.txt
cat header.txt | sed "s/#rname/$refname11/" > ${refname11}_stats.txt
cat header.txt | sed "s/#rname/$refname12/" > ${refname12}_stats.txt
cat header.txt | sed "s/#rname/$refname13/" > ${refname13}_stats.txt
cat header.txt | sed "s/#rname/$refname14/" > ${refname14}_stats.txt
cat header.txt | sed "s/#rname/$refname15/" > ${refname15}_stats.txt

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775
do 
samtools coverage ${sample}_maponly.bam | grep $ref1 | sed "s/${ref1}/$sample/" >> ${refname1}_stats.txt
samtools coverage ${sample}_maponly.bam | grep $ref2 | sed "s/${ref2}/$sample/" >> ${refname2}_stats.txt
samtools coverage ${sample}_maponly.bam | grep $ref3 | sed "s/${ref3}/$sample/" >> ${refname3}_stats.txt
samtools coverage ${sample}_maponly.bam | grep $ref4 | sed "s/${ref4}/$sample/" >> ${refname4}_stats.txt
samtools coverage ${sample}_maponly.bam | grep $ref5 | sed "s/${ref5}/$sample/" >> ${refname5}_stats.txt
samtools coverage ${sample}_maponly.bam | grep $ref6 | sed "s/${ref6}/$sample/" >> ${refname6}_stats.txt
samtools coverage ${sample}_maponly.bam | grep $ref7 | sed "s/${ref7}/$sample/" >> ${refname7}_stats.txt
samtools coverage ${sample}_maponly.bam | grep $ref8 | sed "s/${ref8}/$sample/" >> ${refname8}_stats.txt
samtools coverage ${sample}_maponly.bam | grep $ref9 | sed "s/${ref9}/$sample/" >> ${refname9}_stats.txt
samtools coverage ${sample}_maponly.bam | grep $ref10 | sed "s/${ref10}/$sample/" >> ${refname10}_stats.txt
samtools coverage ${sample}_maponly.bam | grep $ref11 | sed "s/${ref11}/$sample/" >> ${refname11}_stats.txt
samtools coverage ${sample}_maponly.bam | grep $ref12 | sed "s/${ref12}/$sample/" >> ${refname12}_stats.txt
samtools coverage ${sample}_maponly.bam | grep $ref13 | sed "s/${ref13}/$sample/" >> ${refname13}_stats.txt
samtools coverage ${sample}_maponly.bam | grep $ref14 | sed "s/${ref14}/$sample/" >> ${refname14}_stats.txt
samtools coverage ${sample}_maponly.bam | grep $ref15 | sed "s/${ref15}/$sample/" >> ${refname15}_stats.txt
done

cat ${refname1}_stats.txt ${refname2}_stats.txt ${refname3}_stats.txt ${refname4}_stats.txt \
    ${refname5}_stats.txt ${refname6}_stats.txt ${refname7}_stats.txt ${refname8}_stats.txt \
    ${refname9}_stats.txt ${refname10}_stats.txt ${refname11}_stats.txt ${refname12}_stats.txt \
    ${refname13}_stats.txt ${refname14}_stats.txt ${refname15}_stats.txt \
    > all_mapping_stats2.txt
