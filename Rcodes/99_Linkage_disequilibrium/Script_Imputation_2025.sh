
#### This script comes after Main_script and after the Isotype calling in R ####
#### Isotypes_FM scripts create removed_isotypes_Cemee_2025.txt

#mkdir Final_VCF

## Remove the additional lines - I also removed CA2100L45_1_G7 and CA2100L45_2_H7 that were not consistent
bcftools view -S ^removed_isotypes_Cemee_2025.txt Filtered/combined.final.vcf.gz  -O z > Final_VCF/Cemee2025.vcf.gz

### Compute the stats
bcftools stats  Final_VCF/Cemee2025.vcf.gz

### Now


#################### IMPUTATION #########################
#mkdir Final_VCF_output
#mkdir temp_by_chr
## We will imput the lines all at once
#chrom="I"
java -Xmx10g -jar /mnt/data1/mallard/beagle.29Oct24.c8e.jar excludemarkers="Filtered/SNP_2exclude.txt" gt=Final_VCF/Cemee2025.vcf.gz out="temp_by_chr/Cemee2025.imputed_I.vcf.gz" chrom="I" ne=1000 window=5 overlap=2 impute=true imp-segment=0.5 imp-step=0.01 cluster=0.0005 >   /output_impute_chrI.txt &

#chrom="II"
java -Xmx10g -jar /mnt/data1/mallard/beagle.29Oct24.c8e.jar excludemarkers="Filtered/SNP_2exclude.txt" gt=Final_VCF/Cemee2025.vcf.gz out="temp_by_chr/Cemee2025.imputed_II.vcf.gz" chrom="II" ne=1000 window=5 overlap=2 impute=true imp-segment=0.5 imp-step=0.01 cluster=0.0005 > Final_VCF_output/output_impute_chrII.txt &

#chrom="III"
java -Xmx10g -jar /mnt/data1/mallard/beagle.29Oct24.c8e.jar excludemarkers="Filtered/SNP_2exclude.txt" gt=Final_VCF/Cemee2025.vcf.gz out="temp_by_chr/Cemee2025.imputed_III.vcf.gz" chrom="III" ne=1000 window=5 overlap=2 impute=true imp-segment=0.5 imp-step=0.01 cluster=0.0005 > Final_VCF_output/output_impute_chrIII.txt

#chrom="IV"
java -Xmx10g -jar /mnt/data1/mallard/beagle.29Oct24.c8e.jar excludemarkers="Filtered/SNP_2exclude.txt" gt=Final_VCF/Cemee2025.vcf.gz out="temp_by_chr/Cemee2025.imputed_IV.vcf.gz" chrom="IV" ne=1000 window=5 overlap=2 impute=true imp-segment=0.5 imp-step=0.01 cluster=0.0005 > Final_VCF_output/output_impute_chrIV.txt &

#chrom="V"
java -Xmx10g -jar /mnt/data1/mallard/beagle.29Oct24.c8e.jar excludemarkers="Filtered/SNP_2exclude.txt" gt=Final_VCF/Cemee2025.vcf.gz out="temp_by_chr/Cemee2025.imputed_V.vcf.gz" chrom="V" ne=1000 window=5 overlap=2 impute=true imp-segment=0.5 imp-step=0.01 cluster=0.0005 > Final_VCF_output/output_impute_chrV.txt

#chrom="X"
java -Xmx10g -jar /mnt/data1/mallard/beagle.29Oct24.c8e.jar excludemarkers="Filtered/SNP_2exclude.txt" gt=Final_VCF/Cemee2025.vcf.gz out="temp_by_chr/Cemee2025.imputed_X.vcf.gz" chrom="X" ne=1000 window=5 overlap=2 impute=true imp-segment=0.5 imp-step=0.01 cluster=0.0005 > Final_VCF_output/output_impute_chrX.txt &

#chrom="MtDNA"
java -Xmx10g -jar /mnt/data1/mallard/beagle.29Oct24.c8e.jar excludemarkers="Filtered/SNP_2exclude.txt" gt=Final_VCF/Cemee2025.vcf.gz out="temp_by_chr/Cemee2025.imputed_MtDNA.vcf.gz" chrom="MtDNA" ne=1000 window=5 overlap=2 impute=true imp-segment=0.5 imp-step=0.01 cluster=0.0005 > Final_VCF_output/output_impute_chrMtDNA.txt

## Concat all
bcftools concat --naive temp_by_chr/Cemee2025.imputed_*.vcf.gz -O z > Final_VCF/Cemee2025.imputed.vcf.gz
bcftools index Final_VCF/Cemee2025.imputed.vcf.gz

# Name SNPs
path_ref='/projects/nematode_genetics/reference_genomes/20220216_c_elegans_WS283.genome.fa'
bcftools reheader -f ${path_ref}.fai Final_VCF/Cemee2025.imputed.vcf.gz > Final_VCF/Cemee2025.imputed.H.vcf.gz
bcftools annotate --set-id +'%CHROM\_%POS' Final_VCF/Cemee2025.imputed.H.vcf.gz  -O z > Final_VCF/Cemee2025.imputed.H.a.vcf.gz
#bash script_filtration_imputation.sh imputed_All_$pop.ID.vcf.gz $pop

path_gatk='/mnt/data2/desmarais/gatk-4.5.0.0/gatk'

bcftools view -v snps Final_VCF/Cemee2025.imputed.H.a.vcf.gz -O z > Final_VCF/Cemee2025.imputed.SNP.vcf.gz
$path_gatk IndexFeatureFile -I Final_VCF/Cemee2025.imputed.SNP.vcf.gz

# Spot heterozygotes
$path_gatk VariantFiltration -V Final_VCF/Cemee2025.imputed.SNP.vcf.gz --genotype-filter-expression "isHet == 1" --genotype-filter-name "is_het" -O Final_VCF/Cemee2025.imputed.SNP.het.vcf.gz

# Apply high heterozygosity filter
bcftools filter --soft-filter='high_heterozygosity' --mode + --include 'F_PASS(GT="het")<=0.15' -O z Final_VCF/Cemee2025.imputed.SNP.het.vcf.gz > Final_VCF/Cemee2025.imputed.SNP.het1.vcf.gz
$path_gatk IndexFeatureFile -I Final_VCF/Cemee2025.imputed.SNP.het1.vcf.gz

#Put ./. ie missing value for  heterozygotes genotype
$path_gatk SelectVariants \
	-V Final_VCF/Cemee2025.imputed.SNP.het1.vcf.gz \
	--set-filtered-gt-to-nocall \
	-O Final_VCF/Cemee2025.imputed.SNP.without_hetero.vcf.gz

# Apply high missing filter
bcftools filter --soft-filter='high_missing' --mode + --include 'F_MISSING<=0.95' Final_VCF/Cemee2025.imputed.SNP.without_hetero.vcf.gz -O z >  Final_VCF/Cemee2025.imputed.SNP.filtered.vcf.gz

# Hard filtering, remove sites : nonPass, multiallelic, monomorphic
bcftools view -m2 -M2 -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' -f PASS Final_VCF/Cemee2025.imputed.SNP.filtered.vcf.gz -O z > Final_VCF/Cemee2025.imputed.SNP.hard_filtered.vcf.gz

vcftools --maf 0.05 --gzvcf Final_VCF/Cemee2025.imputed.SNP.hard_filtered.vcf.gz --out Final_VCF/Cemee2025.imputed.SNP.filtred.final.vcf --recode
mv Final_VCF/Cemee2025.imputed.SNP.filtred.final.vcf.recode.vcf Final_VCF/Cemee2025.imputed.SNP.filtered.final.vcf
bgzip Final_VCF/Cemee2025.imputed.SNP.filtered.final.vcf


## Rename the lines
bcftools query -l Final_VCF/Cemee2025.imputed.SNP.filtered.final.vcf.gz > All_sample_names.txt
## Change the names in R
#a=read.table("All_sample_names.txt")
#a$V1 <- gsub('CeMee', '', a$V1)
#a$V1 <- gsub('_sorted', '', a$V1)
#write.table(a$V1,file='All_sample_renamed.txt',quote=FALSE,row.names=FALSE,col.names=FALSE)

bcftools reheader -s All_sample_renamed.txt Final_VCF/Cemee2025.imputed.SNP.filtered.final.vcf.gz -o Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz
# Check the results
bcftools query -l Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz > All_sample_names_check.txt

### Archive relevant vcf files ###

cp /mnt/data2/mallard/2025_Cemee_panel_update/Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz /projects/nematode_genetics/Cemee_data/vcf_2025/

### Extract position and chromosome
bcftools query -f '%CHROM\t%POS\n' Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz > Final_VCF/Cemee2025.allpositions.txt

### Stats of the filtered vcf
bcftools stats Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz > Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.stats.txt 

## From Alix github directory
# create chrom map 
mkdir Pruned
bcftools view -H  Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz | cut -f 1 | uniq | awk '{print $0"\t"$0}' > Pruned/chrom-map.txt

### Compute linkage in the Cemee Panel
vcftools --gzvcf /projects/nematode_genetics/Cemee_data/vcf_2025/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --geno-r2 --ld-window-bp 100000 --out ld_decay

#### VCFTOOLS ONLY ON THOR ####
# convert vcf to plink format
vcftools --gzvcf Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --plink --out Pruned/plink.data --chrom-map Pruned/chrom-map.txt
 
#### Linkage disequilibrium ####
vcftools --gzvcf Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --hap-r2 --ld-window-bp 30000 --out Final_VCF/ld_window_30000

# pruning
for r2 in 0.9999 0.99 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
do
	/mnt/data2/desmarais/tool_plink//plink --vcf Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --indep-pairwise 200 10 $r2 --out Pruned/pruned.$r2 --allow-extra-chr --double-id
	/mnt/data2/desmarais/tool_plink/plink --vcf Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --extract Pruned/pruned.$r2.prune.in --recode vcf --out Pruned/pruned.$r2 --allow-extra-chr --double-id # --update-name --set-missing-var-ids @:# 
	bgzip Pruned/pruned.$r2.vcf
done

### Compute linkage blocks with plink ###
mkdir linkage_blocks
/mnt/data2/desmarais/tool_plink/plink --vcf Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --blocks no-pheno-req --maf 0.05 --blocks-max-kb 500 --allow-extra-chr > linkage_blocks/Cemee2025.imputed.SNP.filtered.final_renamed.txt

mkdir linkage_blocks_nomax
cd linkage_blocks_nomax
/mnt/data2/desmarais/tool_plink/plink --vcf ../Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --blocks no-pheno-req --maf 0.05 --allow-extra-chr --blocks-max-kb 20000 > Cemee2025.imputed.SNP.filtered.final_renamed.txt

bcftools view -T sign_positions.txt -O z Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz > Final_VCF/Cemee2025.sign.GWAS.filt.vcf.gz

#mkdir linkage_blocks_nomax_filt
#cd linkage_blocks_nomax_filt
#/mnt/data2/desmarais/tool_plink/plink --vcf ../Final_VCF/Cemee2025.sign.GWAS.filt.vcf.gz --blocks no-pheno-req no-small-max-span --maf 0.05 --allow-extra-chr --blocks-max-kb 500 > Cemee2025.sign.GWAS.filt.txt

#### Linkage disequilibrium ####
mkdir linkage_diseq_filt
cd linkage_diseq_filt
vcftools --gzvcf ../Final_VCF/Cemee2025.sign.GWAS.filt.vcf.gz --hap-r2 --ld-window 100 --out ld_window_100


#/mnt/data2/desmarais/tool_plink/plink --vcf Pruned/pruned.0.9.vcf.gz --blocks no-pheno-req --maf 0.05 --blocks-max-kb 500 --allow-extra-chr > linkage_blocks/pruned.0.9.txt
### Warning: this creates some plink.XX files that are erased progressively

### Compute linkage blocks with plink restricting to the A6140 and CA populations ? ###

###### KINSHIP MATRIX ########
# Rename the chromosome in the vcf file #
bcftools annotate --rename-chrs Pruned/chrom-map_renamed.txt -Oz -o Final_VCF_chr_renamed/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz
## Prune
r2=0.9
	/mnt/data2/desmarais/tool_plink//plink --vcf Final_VCF_chr_renamed/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --indep-pairwise 200 10 $r2 --out Pruned_chr_renamed/pruned.$r2 --allow-extra-chr --double-id
	/mnt/data2/desmarais/tool_plink/plink --vcf Final_VCF_chr_renamed/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --extract Pruned_chr_renamed/pruned.$r2.prune.in --recode vcf --out Pruned_chr_renamed/pruned.$r2 --allow-extra-chr --double-id
	bgzip Pruned_chr_renamed/pruned.$r2.vcf

bcftools query -f '%CHROM\n' Pruned_chr_renamed/pruned.$r2.vcf.gz | sort | uniq

### Make binary on the pruned data for the kinship matrix
mkdir pink_binary_pruned_09
vcftools --gzvcf Pruned_chr_renamed/pruned.0.9.vcf.gz --plink --out pink_binary_pruned_09/plink.data.pruned09 --chrom-map Pruned/chrom-map_renamed.txt
/mnt/data2/desmarais/tool_plink/plink --file pink_binary_pruned_09/plink.data.pruned09 --make-bed --allow-extra-chr --out pink_binary_pruned_09/plink_binary.pruned09

### Create the kinship matrix using GEMMA
mkdir kinship_matrix
## Rename the chromosome for GEMMA
#cut -f1 pink_binary_pruned_09/plink_binary.pruned09.bim | sort | uniq > Pruned/chrom-map_renamed.txt
#/mnt/data2/desmarais/tool_plink/plink --bfile pink_binary_pruned_09/plink_binary.pruned09 --update-chr Pruned/chrom-map_renamed.txt  --make-bed --out pink_binary_pruned_09/plink_binary.pruned09_renamed --allow-extra-chr

## Fix sex info 
awk '{$5=2; $6=1; print}' pink_binary_pruned_09/plink_binary.pruned09.fam > pink_binary_pruned_09/plink_binary.pruned09_fixed.fam
mv pink_binary_pruned_09/plink_binary.pruned09_fixed.fam pink_binary_pruned_09/plink_binary.pruned09.fam
gemma -bfile pink_binary_pruned_09/plink_binary.pruned09 -gk 1 -o kinship_matrix
mv output/kinship_matrix.* kinship_matrix/
/mnt/data2/desmarais/tool_plink/plink --bfile pink_binary_pruned_09/plink_binary.pruned09 --freq

### Make binary plink files for Gemma ###

mkdir pink_binary
/mnt/data2/desmarais/tool_plink/plink --file Pruned/plink.data --make-bed --allow-extra-chr --out pink_binary/plink_binary

### Run Gemma
###
# Use the latest installed gemma version
conda activate gemma-env 

### Only takes one column at a time
cut -f 1 phenotypic_data/gemma_pheno.txt > phenotypic_data/gemma_pheno_T1.txt
nohup gemma -bfile pink_binary/plink_binary -k kinship_matrix/kinship_matrix.cXX.txt -lmm 4 -p phenotypic_data/gemma_pheno_T1.txt -o gwas_results_T1 &

for i in 8 9 10 11 12 13 14
do
cut -f $i phenotypic_data/gemma_pheno.txt > phenotypic_data/gemma_pheno_T${i}.txt
nohup gemma -bfile pink_binary/plink_binary -k kinship_matrix/kinship_matrix.cXX.txt -lmm 4 -p phenotypic_data/gemma_pheno_T${i}.txt -o gwas_results_T${i} > gemma_output_txt/output_gemma_t${i}.txt &
done

######## Could we do permutations ? How long would it take ?

nohup gemma -bfile pink_binary/plink_binary -k kinship_matrix/kinship_matrix.cXX.txt -lmm 1 -p phenotypic_data/gemma_pheno_T1.txt -o gwas_results_T1 &


##### For Stephanie/Itamar (and probably all PoolSeq data analysis): keep all Cemee variants, do not filter for maf 5%
vcftools --maf 0.00001 --gzvcf Final_VCF/Cemee2025.imputed.SNP.hard_filtered.vcf.gz --out Final_VCF/Cemee2025.imputed.SNP.filtred.final_noMAF_filtering.vcf --recode
mv Final_VCF/Cemee2025.imputed.SNP.filtred.final_noMAF_filtering.vcf.recode.vcf Final_VCF/Cemee2025.imputed.SNP.filtered.final_noMAF_filtering.vcf
bgzip Final_VCF/Cemee2025.imputed.SNP.filtered.final_noMAF_filtering.vcf


## Rename the lines ## Should all be the same
bcftools reheader -s All_sample_renamed.txt Final_VCF/Cemee2025.imputed.SNP.filtered.final_noMAF_filtering.vcf.gz -o Final_VCF/Cemee2025.imputed.SNP.filtered.final_noMAF_filtering_renamed.vcf.gz

### Archive relevant vcf files ###
cp /mnt/data2/mallard/2025_Cemee_panel_update/Final_VCF/Cemee2025.imputed.SNP.filtered.final_noMAF_filtering_renamed.vcf.gz /projects/nematode_genetics/Cemee_data/vcf_2025/

### Extract position and chromosome
bcftools query -f '%CHROM\t%POS\n' Final_VCF/Cemee2025.imputed.SNP.filtered.final_noMAF_filtering_renamed.vcf.gz > Final_VCF/Cemee2025.allpositions_noMAF.txt

### Stats of the filtered vcf
bcftools stats Final_VCF/Cemee2025.imputed.SNP.filtered.final_noMAF_filtering_renamed.vcf.gz > Final_VCF/Cemee2025.imputed.SNP.filtered.final_renamed_noMAF.stats.txt 
 

bcftools stats Final_VCF/Cemee2025.vcf.gz > test_FM.txt


