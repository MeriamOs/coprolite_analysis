#!/bin/bash

input_tped="angsd_bam_trimmed.tped"
output_list="multiallelic_sites.txt"

awk '
{
    n = split("", alleles)  # clear
    delete alleles
    for (i = 5; i <= NF; i++) {
        if ($i != "0") alleles[$i]++
    }
    if (length(alleles) > 2) {
        print $2  # print SNP ID (2nd column)
    }
}' "$input_tped" > "$output_list"

echo "Done. Multiallelic site list saved to $output_list"
