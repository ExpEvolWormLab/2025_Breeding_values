

## Try to load the gff3

library(rtracklayer)
library(GenomicRanges)

#gff <- import("/Volumes/nematode_genetics/reference_genomes/c_elegans.PRJNA13758.WS283.annotations.gff3.gz")
gff <- import("/Volumes/nematode_genetics/reference_genomes/WS283_gene_only_col2.gff2")
str(gff)
genes <- gff[gff$type == "gene"]
head(genes)


#all_cluster_merged_with_selected_groups = subset(all_cluster_merged,group%in%subset(all_GWAS_clusters,is_sign_T7)$group)
#dim(all_cluster_merged_with_selected_groups) # 5361

#df=subset(all_cluster_merged_with_selected_groups,is_sign_T7)[,c("CHROM","POS")]

df=subset(all_cluster_merged,is_sign_T7)[,c("CHROM","POS")]
df$CHROMOSOME=paste0("CHROMOSOME_",df$CHROM)

gr_pos <- GRanges(seqnames = df$CHROMOSOME,
                  ranges = IRanges(df$POS, df$POS))
hits <- findOverlaps(gr_pos, genes)

result <- data.frame(
  CHROM = df$CHROMOSOME[queryHits(hits)],
  POS   = df$POS[queryHits(hits)],
  Gene  = genes$Gene[subjectHits(hits)]
)

dim(result)
length(unique(result$Gene))
write.table(unique(result$Gene),file="test.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

### Intercept T7/T14 ?

#df=subset(all_cluster_merged_with_selected_groups,is_sign_T7 & is_sign_T14)[,c("CHROM","POS")]
df=subset(all_cluster_merged,is_sign_T7 & is_sign_T14)[,c("CHROM","POS")]
df$CHROMOSOME=paste0("CHROMOSOME_",df$CHROM)

gr_pos <- GRanges(seqnames = df$CHROMOSOME,
                  ranges = IRanges(df$POS, df$POS))
hits <- findOverlaps(gr_pos, genes)

result <- data.frame(
  CHROM = df$CHROMOSOME[queryHits(hits)],
  POS   = df$POS[queryHits(hits)],
  Gene  = genes$Gene[subjectHits(hits)]
)

dim(result)
length(unique(result$Gene))
write.table(unique(result$Gene),file="test.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)



load('RData/Genotype_info_Cemee.RData')
snp_positions=as.data.frame(snp_positions)
head(snp_positions)


df=snp_positions
df$CHROMOSOME=paste0("CHROMOSOME_",df$CHROM)
df$POS=as.numeric(df$POS)
gr_pos <- GRanges(seqnames = df$CHROMOSOME,
                  ranges = IRanges(df$POS, df$POS))
hits <- findOverlaps(gr_pos, genes)

result <- data.frame(
  CHROM = df$CHROMOSOME[queryHits(hits)],
  POS   = df$POS[queryHits(hits)],
  Gene  = genes$Gene[subjectHits(hits)]
)

dim(result)
length(unique(result$Gene))
write.table(unique(result$Gene),file="bkgd.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)


