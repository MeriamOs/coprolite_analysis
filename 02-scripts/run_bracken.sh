module purge
module load Bracken/2.7-GCC-11.3.0
module load Python/3.11.3-gimkl-2022a

#ERR5863536 SRR12455959 ERR3761407 ERR3761412 SRR14842400 SRR14842401 SRR14842402 SRR14842403 SRR14842404 SRR14842405 SRR14842406 SRR14842408 SRR14842409 SRR14842410 SRR14842411 SRR14842412 SRR14842413 SRR14842414 SRR14842416 SRR14842417 SRR14842418 SRR14842419 SRR14842320 SRR14842420 SRR14842321 SRR14842421 SRR14842322 SRR14842422 SRR14842323 SRR14842423 SRR14842425 SRR14842326 SRR14842329 SRR14842330 SRR14842331 SRR14842332 SRR14842333 SRR14842339 SRR14842343 SRR14842349 SRR14842353 SRR14842354 SRR14842355 SRR14842356 SRR14842357 SRR14842358 SRR14842359 SRR14842360 SRR14842361 SRR14842362 SRR14842364 SRR14842365 SRR14842366 SRR14842367 SRR14842368 SRR14842369 SRR14842370 SRR14842371 SRR14842372 SRR14842373 SRR14842375 SRR14842376 SRR14842377 SRR14842378 SRR14842379 SRR14842380 SRR14842381 SRR14842382 SRR14842383 SRR14842384 SRR14842386 SRR14842387 SRR14842388 SRR14842389 SRR14842390 SRR14842391 SRR14842392 SRR14842393 SRR14842394 SRR14842395 SRR14842397 SRR14842398 SRR14842399 SRR14842324 SRR14842325 SRR14842327 SRR14842328 SRR14842334 SRR14842335 SRR14842336 SRR14842337 SRR14842338 SRR14842340 SRR14842341 SRR14842342 SRR14842344 SRR14842345 SRR14842346 SRR14842347 SRR14842348 SRR14842350 SRR14842351 SRR14842352 SRR14842363 SRR14842374 SRR14842385 SRR14842396 SRR14842407 SRR14842415 SRR14842424 SRR14842426
#MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775

DB='/nesi/nobackup/uoo02328/meriam/kraken/databases/2024-standard'
#'/nesi/nobackup/uoo02328/meriam/kraken/databases/2024-pluspfp'
prefix=2024-standard
#2024-pluspfp

for sample in MS10790 MS10902 MS10903 MS10904 MS11102 MS11103 MS11107 Blank1_WH Blank2_WH KH_blank_1 KH_blank_2 LB_blank_1 LB_blank_2;
do
bracken -d $DB -i reports/${sample}_kraken_${prefix}.50.txt -o bracken_phylum/${sample}.bracken -l P
python ../KrakenTools/filter_bracken.out.py -i bracken_phylum/${sample}.bracken --exclude 7711 -o bracken_phylum/${sample}_filtered.bracken
done

# TAXIDs
# Species Homo sapiens 9606
# Genus Homo 9605
# Phylum Chordata 7711 

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
do
bracken -d $DB -i reports_modern/${sample}_kraken_${prefix}_conf0.50.txt -o bracken_genus/${sample}.bracken -l G
python ../KrakenTools/filter_bracken.out.py -i bracken_genus/${sample}.bracken --exclude 9605 -o bracken_genus/${sample}_filtered.bracken
done

#Blank1_WH Blank2_WH KH_blank_1 KH_blank_2 LB_blank_1 LB_blank_2
for sample in MS10790 MS10902 MS10903 MS10904 MS11102 MS11103 MS11107 MS11108 Blank1_WH Blank2_WH KH_blank_1 KH_blank_2 LB_blank_1 LB_blank_2;
do
bracken -d $DB -i reports/${sample}_kraken_${prefix}_conf0.50.txt -o bracken_decontam/${sample}.bracken -l S
python ../KrakenTools/filter_bracken.out.py -i bracken_decontam/${sample}.bracken --exclude 9606 -o bracken_decontam/${sample}_filtered.bracken
done

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
do
bracken -d $DB -i reports/${sample}_kraken_${prefix}_conf0.50.txt -o bracken_species/${sample}.bracken -l S
python ../KrakenTools/filter_bracken.out.py -i bracken_species/${sample}.bracken --exclude 9606 -o bracken_species/${sample}_filtered.bracken
done

python combine_bracken_outputs.py \
--files bracken_phylum/*_filtered.bracken \
-o palaeofaeces_modern_combined_bracken_phylum.bracken

python combine_bracken_outputs.py \
--files bracken_genus/*_filtered.bracken \
-o palaeofaeces_modern_combined_bracken_genus.bracken

python combine_bracken_outputs.py \
--files bracken_decontam/*_filtered.bracken \
-o negatives_combined_${prefix}_bracken_species.bracken

python combine_bracken_outputs.py \
--files bracken_species/*_filtered.bracken \
-o palaeofaeces_combined_${prefix}_bracken_species.bracken

#python tsv_to_csv.py < input_file.tsv > output_file.csv
#python tsv_to_csv.py < palaeofaeces_modern_combined_bracken_phylum.bracken > palaeofaeces_modern_combined_bracken_phylum.csv
#python tsv_to_csv.py < palaeofaeces_modern_combined_bracken_genus.bracken > palaeofaeces_modern_combined_bracken_genus.csv

python tsv_to_csv.py < palaeofaeces_combined_${prefix}_bracken_species.bracken > palaeofaeces_combined_${prefix}_bracken_species.csv
python tsv_to_csv.py < negatives_combined_${prefix}_bracken_species.bracken > negatives_combined_${prefix}_bracken_species.csv