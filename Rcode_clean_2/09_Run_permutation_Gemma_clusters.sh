#!/bin/bash

for j in 2 4 6 8 9 10
do
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
do
cut -f $i data_BLUE/gemma_pheno.txt > data_BLUE/split/gemma_pheno_T${i}.txt
gemma -g data_BLUE_clusters/genotypes_clusters_R${j}.txt -k data_BLUE/kinship.txt -lmm 4 -p data_BLUE/split/gemma_pheno_T${i}.txt -o gwas_results_clusters_BLUE_R${j}_T${i} > gemma_output_BLUE_clusters_txt/output_gemma_clusters_R${j}_t${i}.txt
done
wait
done
