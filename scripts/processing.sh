#!/bin/bash

# Identify nitrogen cycle genes across all metagenomes: Subset gene phylogeny and product abundance
for dir in *metagenome*;do
    cd $dir
    cd IMG\ Data/
    bash ../../script2.sh
    cd ../../
done

# Attach metagenome name to gene phylogeny and product abundance files
for dir in *metagenome*;do
    cd $dir
    for i in $(ls *_nitrogen_gene_phylogeny.txt);do
        mv $i ${i%_nitrogen_gene_phylogeny*}_nitrogen_gene_phylogeny_$dir.txt
    done
    for i in $(ls *_nitrogen_product_abundance.txt);do
        mv $i ${i%_nitrogen_product_abundance*}_nitrogen_product_abundance_$dir.txt
    done
    cd ../
done

# Move each phylogeny and product abundance file into the main directory
for dir in *metagenome*;do
    cd $dir
    for i in $(ls *_nitrogen_gene_phylogeny*);do
        mv $i ../
    done
    for i in $(ls *_nitrogen_product_abundance*);do
        mv $i ../ 
    done
    cd ../
done
