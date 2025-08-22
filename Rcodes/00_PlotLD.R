
### Broken code, needs fix

rm(list=ls())
gc()

library("ggplot2")
ld <- read.table("ld_window_30000.hap.ld", header=TRUE,nrows = 10000)
ld$distance <- ld$POS2 - ld$POS1
summary((ld$distance))
ggplot(ld, aes(x = distance, y = R.2)) +
  geom_point(alpha = 0.3, size = 1) +
  labs(x = "Distance between SNPs (bp)", y = expression(r^2),
       title = "LD Decay Plot") +
  theme_minimal()+geom_smooth(method = "loess")

library(SNPRelate)
# Convert VCF to GDS (fast binary storage)
#snpgdsVCF2GDS("/Volumes/nematode_genetics/Cemee_data/vcf_2025/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz", "data.gds", method="biallelic.only")

# Open the GDS file
genofile <- snpgdsOpen("data.gds")

# SNP annotations
snp_annot <- snpgdsSNPList(genofile) %>%
  select(snp.id, chromosome, position)

snp.id.CHRV <- subset(snp_annot,chromosome=='V')$snp.id
  
ld <- snpgdsLDMat(genofile, method="r", slide=100, 
                  #ld.threshold=NA, 
                  sample.id=NULL, 
                  #maf=0.05, 
                  #missing.rate=0.05, 
                  num.thread=4, 
                  verbose=TRUE,
                  snp.id=snp.id.CHRV)
# Close GDS file
snpgdsClose(genofile)

# Extract SNP IDs and positions
snp_ids <- ld$snp.id
snp_pos <- snp_annot$position[match(snp_ids, snp_annot$snp.id)]

# Flatten upper triangle of LD matrix into a tidy data frame
ld_df <- which(upper.tri(ld$LD), arr.ind = TRUE) %>%
  as.data.frame() %>%
  mutate(
    SNP1 = snp_ids[row],
    SNP2 = snp_ids[col],
    Pos1 = snp_pos[row],
    Pos2 = snp_pos[col],
    Dist = abs(Pos1 - Pos2),
    R2 = (ld$LD)[cbind(row, col)]^2
  )


# Prepare data frame for plotting
#library(dplyr)
#ld_df <- data.frame(SNP1 = ld$snp.id[ld$snp.id1],
#                    SNP2 = ld$snp.id[ld$snp.id2],
#                    R2 = ld$ld^2) %>%
#  left_join(snp_annot, by=c("SNP1"="snp.id")) %>%
#  left_join(snp_annot, by=c("SNP2"="snp.id"))

dim(ld_df)

# Plot LD decay
#library(ggplot2)
ggplot(ld_df[sample(1:nrow(ld_df),10000),], aes(x=abs(Pos1 - Pos2), y=R2)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_smooth() +
  labs(x="Distance (bp)", y=expression(r^2), title="LD decay") +
  theme_minimal()
