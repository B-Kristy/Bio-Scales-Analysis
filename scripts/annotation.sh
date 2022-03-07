#!/bin/bash

# Remove spaces in phylogeny and product files
for i in $(ls *_product_names.tsv);do
     awk '{gsub(" ",";",$0); print}' $i | awk -F '\t' '{print $1"\t"$2}' > ${i%_product_names*}_product_names2.tsv
done

for i in $(ls *_gene_phylogeny.tsv);do
    awk '{gsub(" ",";",$0); print}' $i | awk -F '\t' '{print $1"\t"$5}' > ${i%_gene_phylogeny*}_gene_phylogeny2.tsv
done

# Append gene phylogeny with product names
for i in $(ls *_product_names2.tsv);do
    awk -F '\t' 'FNR==NR{a[$1]=$2 FS $3;next}{ print $0, a[$1]}' OFS='\t' $i ${i%_product_names2*}_gene_phylogeny2.tsv > ${i%_product_names2*}_test.txt
done

# Append contig names to scaffold read depth
for i in $(ls *_contig_names_mapping.tsv);do
    awk -F '\t' 'FNR==NR{a[$1]=$2 FS $3; next}{ print $0, a[$1]}' $i pairedMapped_sorted.bam.cov > ${i%_contig_names_mapping*}_test2.txt
done

# Print contig and read depth and replace space with an underscore
for i in $(ls *_test2.txt);do
    awk 'NR>1{print $12, $2}' $i | awk '{gsub (" ", "_", $0); print}' > ${i%_test2*}_test3.txt
done

# Replace "_" with "\t"
for i in $(ls *_test3.txt);do
    awk '{gsub("_","\t",$0); print}' $i > ${i%_test3*}_test4.txt
done

for i in $(ls *_test.txt);do
    awk '{gsub("_","\t",$0); print}' $i > ${i%_test*}_test5.txt     
done

# Append read coverage to each contig_gene
for i in $(ls *_test4.txt);do
    awk -F '\t' 'FNR==NR{a[$2] = $2 FS $3; next}{ print $0, a[$2]}' OFS='\t' $i ${i%_test4*}_test5.txt > ${i%_test4*}_test6.txt
done

# Filter file based on product name (Subset nitrogen cycle genes)
for i in $(ls *_test6.txt);do
    grep 'nitrogenase\|ammonia;monooxygenase\|hydoxylamine;oxidoreductase\|nitrite;oxidoreductase\|nitrate;reductase\|nitrite;reductase\|nitrous;oxide;reductase\|nitric;oxide;reductase' $i > ${i%_test6*}_test7.txt 
done     

# Split into two files: product taxa and product abundance 
for i in $(ls *_test7.txt); do
    awk -F '\t' '{print $1"_"$2"_"3"_"$4"\t"$5"\t"$6}' $i > ${i%_test7*}_nitrogen_gene_phylogeny.txt
    awk -F '\t' '{print $1"_"$2"_"$3"_"$4"\t"$9}' $i > ${i%_test7*}_nitrogen_product_abundance.txt
done

# Remove all intermediate files
rm *test*
rm *2*

# move output files into main directory
for i in $(ls *nitrogen*);do
    mv $i ../
done