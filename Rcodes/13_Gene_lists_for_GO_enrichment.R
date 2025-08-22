rm(list=ls())
library(vcfR)
library(dplyr)
library(readr)
library(data.table)
library(stringr)
library(data.table)
library(gplots)
library(limma) 

# Load GWAS results
load(file="RData/GWAS_results_concatenated.RData") 
# Load GWAS results for clusters 
load(file="RData/GWAS_cluster_results_concatenated.RData") 

### Sign in NGM
sign_NGM <- rowSums(all_GWAS_clusters[,paste0("is_sign_T",1:7)])>0 # 4171 clus
sign_NaCl <-rowSums(all_GWAS_clusters[,paste0("is_sign_T",8:14)])>0 # 3345 clust

chisq.test(table(sign_NGM,sign_NaCl))$p.value
### Significant dependence
### Trait by trait ??
ind_per_trait <- NULL
for(i in 1:7){
  ind_per_trait=c(ind_per_trait,
                  chisq.test(table(all_GWAS_clusters[,paste0("is_sign_T",i)],all_GWAS_clusters[,paste0("is_sign_T",(i+7))]))$p.value)
}
print(ind_per_trait) # all 0

vennDiagram(vennCounts(
  cbind(sign_NGM,sign_NaCl)),
  names=c("all 7 NGM","all 7 NaCl"),main='All')

#### Export the list of complete list of SNPs
SNP_Cemee=data.frame(all_GWAS[,c(1,3)])
names(SNP_Cemee)=c("CHROM","POS")

SNP_Cemee$CHROM =  (tstrsplit(all_GWAS$rs,'_')[[1]])

#get genes from positions with same code as before  #CHROM & POS (instead of CHROM &BP)
# Load gene list: 
gene_locations <- read.table(
  "data/WS283_gene_only_col2.gff2",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  quote = "",              # prevents quotes from messing with parsing
  stringsAsFactors = FALSE,
  fill = TRUE              # fills missing fields so line length mismatches don't crash
)

#rename columns #chrom start and end is important, the other things chatGPT suggested, don't trust!! 
colnames(gene_locations) <- c(
  "chrom", "source", "type", "start", "end",
  "score", "strand", "phase", "attributes"
)

# Convert to data.tables explicitly
setDT(gene_locations)
#make sure chromosome naming is the same
gene_locations[, CHROM := gsub("^CHROMOSOME_", "", chrom)]
# Prepare genes intervals
genes_intervals <- gene_locations[, .(CHROM, start, end, attributes)]
setkey(genes_intervals, CHROM, start, end)

extract_genes_from_SNPs <- function(SNP_df){
  # Convert to data.tables explicitly
  setDT(SNP_df)
  # Add start/end columns as the SNP position (point intervals)
  SNP_df[, start := POS]
  SNP_df[, end := POS]
  # foverlaps expects keys on the tables
  setkey(SNP_df, CHROM, start, end)

    All_SNP_gene_map <- foverlaps(
    x = SNP_df,
    y = genes_intervals,
    by.x = c("CHROM", "start", "end"),
    by.y = c("CHROM", "start", "end"),
    type = "within",
    nomatch = NULL
  )
  gene_snp_summary <- All_SNP_gene_map[
    , .(
      gene_start   = unique(start)[1],   # gene start from gene_locations
      gene_end     = unique(end)[1],     # gene end from gene_locations
      SNP_count    = .N,
      SNP_positions = paste(unique(i.start), collapse = ", ")  # SNP positions
    ),
    by = attributes  # group by gene
  ]
  
  ###################################
  #now split the content of the attributes column
  #now split the attributes to extra column
  # Extract key-value pairs from attributes
  attr_dt <- tstrsplit(gsub('"', '', gene_snp_summary$attributes), " ; ", fixed = TRUE)
  
  # Split each "key value" into separate named columns
  attr_list <- lapply(attr_dt, function(col) {
    kv <- tstrsplit(col, " ", fixed = TRUE)
    data.table(Key = kv[[1]], Value = kv[[2]])
  })
  
  # Combine all into one data.table
  attr_wide <- rbindlist(attr_list, use.names = TRUE, fill = TRUE)
  attr_wide <- dcast(attr_wide, rowid(Key) ~ Key, value.var = "Value")
  
  # Remove helper column
  attr_wide[, rowid := NULL]
  
  # Bind to original data (dropping old attributes column)
  gene_snp_summary_clean <- cbind(attr_wide, gene_snp_summary)
  return(gene_snp_summary_clean)
  }

#save the background list of genes
all_gene_list <- extract_genes_from_SNPs(SNP_Cemee)
write.table(all_gene_list$Gene, "table_for_GO/All_Cemee_gene.txt", row.names = FALSE, quote = FALSE,col.names=FALSE)
length(unique(all_gene_list$Gene)) # 19k genes

all_GWAS$CHROM =  (tstrsplit(all_GWAS$rs,'_')[[1]])
all_GWAS$POS = as.numeric(tstrsplit(all_GWAS$rs,"_")[[2]])

NGM_all_Traits <- extract_genes_from_SNPs(subset(all_GWAS,rowSums(all_GWAS[,paste0("is_sign_T",1:7)])>0)[,c("CHROM","POS")])
write.table(NGM_all_Traits$Gene, "table_for_GO/NGM_all_traits_gene.txt", row.names = FALSE, quote = FALSE,col.names=FALSE)
length(unique(NGM_all_Traits$Gene)) # 8662

NaCl_all_Traits <- extract_genes_from_SNPs(subset(all_GWAS,rowSums(all_GWAS[,paste0("is_sign_T",8:14)])>0)[,c("CHROM","POS")])
write.table(NaCl_all_Traits$Gene, "table_for_GO/NaCl_all_traits_gene.txt", row.names = FALSE, quote = FALSE,col.names=FALSE)
length(unique(NaCl_all_Traits$Gene)) # 7578

Shared_env_all_Traits <- NaCl_all_Traits$Gene[NaCl_all_Traits$Gene%in%NGM_all_Traits$Gene]
length(Shared_env_all_Traits) # 4895
write.table(Shared_env_all_Traits, "table_for_GO/Shared_all_traits_gene.txt", row.names = FALSE, quote = FALSE,col.names=FALSE)


### Body Area NGM / NaCl
NGM_area <- extract_genes_from_SNPs(subset(all_GWAS,(all_GWAS[,paste0("is_sign_T",7)])>0)[,c("CHROM","POS")])
write.table(NGM_area$Gene, "table_for_GO/NGM_area_gene.txt", row.names = FALSE, quote = FALSE,col.names=FALSE)
length(unique(NGM_area$Gene)) # 2815
NaCl_area <- extract_genes_from_SNPs(subset(all_GWAS,(all_GWAS[,paste0("is_sign_T",14)])>0)[,c("CHROM","POS")])
write.table(NaCl_area$Gene, "table_for_GO/NaCl_area_gene.txt", row.names = FALSE, quote = FALSE,col.names=FALSE)
length(unique(NaCl_area$Gene)) # 2916

### Area shared ?
Shared_area <- NGM_area$Gene[NGM_area$Gene%in%NaCl_area$Gene]
length(Shared_area) # 1209
write.table(Shared_area, "table_for_GO/Share_area_gene.txt", row.names = FALSE, quote = FALSE,col.names=FALSE)


