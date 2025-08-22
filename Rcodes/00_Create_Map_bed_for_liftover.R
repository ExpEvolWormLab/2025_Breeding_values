rm(list=ls())
gc()
library(data.table)
## Load the map to define arm/center/tip

rk_map <- read.table("Caendr_filtering/rockman_kruglya_map.txt",h=TRUE,sep='\t')
head(rk_map)
map <- read.table("Caendr_filtering/c_elegans_genetic_map.tsv",h=TRUE)
### Test the map -> OK
 par(mfrow=c(2,3))
 for(i in 1:6){
 mapI=subset(map,chrom==rk_map$CHROM[i])
 plot(mapI$marker_pos,mapI$genetic,type='l')
 abline(v=cumsum(1000*as.numeric(rk_map[i,2:6])))
 }
i=1
cumsum(1000*as.numeric(rk_map[i,2:6]))
 start_end_coordinates <- apply(as.matrix(rk_map[,2:6]),1,function(x){
   cumsum(1000*as.numeric(x))})
 start_end_coordinates=t(start_end_coordinates)
 
bed_coordinates <- NULL
for(i in 1:nrow(rk_map)){
  for(j in 1:ncol(start_end_coordinates)){
    if(j==1){
    bed_coordinates=c(bed_coordinates,paste0(rk_map$CHROM[i],":",0,"-",start_end_coordinates[i,j]))
    }else{
      bed_coordinates=c(bed_coordinates,paste0("chr",rk_map$CHROM[i],":",start_end_coordinates[i,(j-1)],"-",start_end_coordinates[i,j]))
    }
    
  }
  
}
 data.frame(bed_coordinates)

 #### We then translate from ce6 -> ce11 using :
 #https://genome.ucsc.edu/cgi-bin/hgLiftOver
 rk_map_lift_bed <- read.table("Caendr_filtering/rockman_kruglya_map_liftOver.bed",sep='\t')
 rk_map_lift_bed_df=data.frame(
   CHROM=substring(tstrsplit(rk_map_lift_bed$V1,":")[[1]],4,8),
   Start=c(0,as.numeric(tstrsplit(rk_map_lift_bed$V1,"-")[[2]])[1:29]),
   End=as.numeric(tstrsplit(rk_map_lift_bed$V1,"-")[[2]]),
   type=rep(names(rk_map[2:6]),6))
 rk_map_lift_bed_df$Start[rk_map_lift_bed_df$type=='left_tip']=0
 ## Very little difference :)
 write.table(rk_map_lift_bed_df,file='Caendr_filtering/rockman_kruglya_map_liftOver.txt',sep='\t',row.names=FALSE,quote=FALSE)
 
 
 par(mfrow=c(2,3))
 for(i in 1:6){
   mapI=subset(map,chrom==rk_map$CHROM[i])
   plot(mapI$marker_pos,mapI$genetic,type='l')
   abline(v=cumsum(1000*as.numeric(rk_map[i,2:6])))
   abline(v=subset(rk_map_lift_bed_df,CHROM==rk_map$CHROM[i])$Coord,col="red")
 }
 
 # 
