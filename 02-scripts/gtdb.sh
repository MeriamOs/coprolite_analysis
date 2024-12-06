module purge

module load MetaPhlAn/4.1.0-gimkl-2022a-Python-3.10.5

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 MS11679 MS11683 MS11686 MS11770 MS11771 MS11774 MS11775;
do
sgb_to_gtdb_profile.py -i ${sample}_4.1.0.metaphlan_profile.txt -o ${sample}_gtdb_profile.txt
done

merge_metaphlan_tables.py *gtdb_profile.txt --gtdb_profiles > merged_gtdb_abundance_table.txt

#create species only table
grep -E "s__|MS" merged_gtdb_abundance_table.txt \
| sed "s/^.*;//g" \
| sed "s/MS[0-9]*-//g" \
> merged_gtdb_abundance_table_species.txt

module purge

module load Python/2.7.18-gimkl-2020a

python /home/vanme090/.local/bin/hclust2.py \
-i merged_gtdb_abundance_table_species.txt \
-o metaphlan4_gtdb_abundance_heatmap_species.png \
--f_dist_f braycurtis \
--s_dist_f braycurtis \
--cell_aspect_ratio 0.5 \
--log_scale \
--flabel_size 10 --slabel_size 10 \
--max_flabel_len 100 --max_slabel_len 100 \
--minv 0.1 \
--dpi 300 \
--ftop 25