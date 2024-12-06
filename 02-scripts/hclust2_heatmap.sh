#module load Python/2.7.18-gimkl-2020a
#pip install hclust2


python /home/vanme090/.local/bin/hclust2.py \
-i merged_abundance_table_species.txt \
-o metaphlan4_abundance_heatmap_species.png \
--f_dist_f braycurtis \
--s_dist_f braycurtis \
--cell_aspect_ratio 0.5 \
--log_scale \
--flabel_size 10 --slabel_size 10 \
--max_flabel_len 100 --max_slabel_len 100 \
--minv 0.1 \
--dpi 300 \
--ftop 25
