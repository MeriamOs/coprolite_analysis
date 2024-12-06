module load Python/2.7.18-gimkl-2020a

#pip install graphlan

#pip install export2graphlan

#tail -n +2 merged_abundance_table.txt | cut -f1,3- > merged_abundance_table_reformatted.txt
 
#python /home/vanme090/.local/bin/export2graphlan.py --skip_rows 1 -i merged_abundance_table_reformatted.txt --tree merged_abundance.tree.txt --annotation merged_abundance.annot.txt --most_abundant 100 --abundance_threshold 1 --least_biomarkers 10 --annotations 5,6 --external_annotations 7 --min_clade_size 1

#python /home/vanme090/.local/bin/graphlan_annotate.py --annot merged_abundance.annot.txt merged_abundance.tree.txt merged_abundance.xml

#python /home/vanme090/.local/bin/graphlan.py --dpi 300 merged_abundance.xml merged_abundance.png --external_legends

#tail -n +2 merged_gtdb_abundance_table.txt | cut -f1,3- > merged_gtdb_abundance_table_reformatted.txt
 
python /home/vanme090/.local/bin/export2graphlan.py --skip_rows 1 -i merged_gtdb_abundance_table_reformatted.txt --tree merged_gtdb_abundance.tree.txt --annotation merged_gtdb_abundance.annot.txt --most_abundant 100 --abundance_threshold 1 --least_biomarkers 10 --annotations 5,6 --external_annotations 7 --min_clade_size 1

python /home/vanme090/.local/bin/graphlan_annotate.py --annot merged_gtdb_abundance.annot.txt merged_gtdb_abundance.tree.txt merged_gtdb_abundance.xml

python /home/vanme090/.local/bin/graphlan.py --dpi 300 merged_gtdb_abundance.xml merged_gtdb_abundance.png --external_legends