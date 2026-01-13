rm(list=ls())
library(vcfR)

map <- read.table("map/c_elegans_genetic_map.tsv",sep='\t',h=TRUE)

chr_temp = 'I'

for(chr_temp in c("I","II","III","IV","V","X")){

vcf_I <- read.vcfR(paste0("vcf/A6140_derived.",chr_temp,".vcf.gz"))
vcf_I.fix <- data.frame(vcf_I@fix)
map.I = subset(map,chrom==chr_temp)

## Compute recombination rates in the region 

map.I$rates <- 0
map.I$ends  <- 0

map.I$starts[1]=0
map.I$ends[1]= map.I$marker_pos[1]
map.I$ends[1]= map.I$marker_pos[1]

dcm <- map.I$genetic[2] - map.I$genetic[1]
dbp <- map.I$marker_pos[2] - map.I$marker_pos[1]

map.I$rates[1] <-  (dcm * 0.01) / dbp

for (i in 2:(nrow(map.I))) {
  dcm <- map.I$genetic[i] - map.I$genetic[i-1]
  dbp <- map.I$marker_pos[i] - map.I$marker_pos[i-1]
  map.I$rates[i] <-  (dcm * 0.01) / dbp
  map.I$ends[i]  <- map.I$marker_pos[i]
  map.I$starts[i] <- map.I$marker_pos[i-1] 
}


recom.map <- data.frame(rates=0,ends=vcf_I.fix$POS)
vcf_I.fix$POS=as.numeric(vcf_I.fix$POS)
for(i in 1:nrow(vcf_I.fix)){
  ## Find the closest ends ?
  dist_maps.ends <- abs(vcf_I.fix$POS[i]-map.I$ends)
  dist_maps.starts <- abs(vcf_I.fix$POS[i]-map.I$ends)
if(min(dist_maps.ends)<min(dist_maps.starts)){
  recom.map$rates[i]=map.I$rates[dist_maps.ends==min(dist_maps.ends)]
}else{
  recom.map$rates[i]=map.I$rates[dist_maps.starts==min(dist_maps.starts)]
}
}

## Some issues with extreme values ?
options(scipen = 999)

write.table(recom.map[,1],file=paste0('map/CHR',chr_temp,'.map.rates.txt'),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(recom.map[,2],file=paste0('map/CHR',chr_temp,'.map.ends.txt'),quote=FALSE,row.names=FALSE,col.names=FALSE)

}

### Remove missing data:
#bcftools view -i 'F_MISSING=0' A6140_derived.I.vcf.gz -o A6140_derived.I_nomissing.vcf
#bcftools view -i 'F_MISSING=0' A6140_derived.II.vcf.gz -o A6140_derived.II_nomissing.vcf
#bcftools view -i 'F_MISSING=0' A6140_derived.III.vcf.gz -o A6140_derived.III_nomissing.vcf
#bcftools view -i 'F_MISSING=0' A6140_derived.IV.vcf.gz -o A6140_derived.IV_nomissing.vcf
#bcftools view -i 'F_MISSING=0' A6140_derived.V.vcf.gz -o A6140_derived.V_nomissing.vcf
#bcftools view -i 'F_MISSING=0' A6140_derived.X.vcf.gz -o A6140_derived.X_nomissing.vcf

