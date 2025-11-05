

rm(list=ls());gc()
library(ggplot2)
library(dplyr)

ld=read.table("/Users/mallard/Documents/GitHub/2025_Breeding_values/Rcodes/99_Linkage_disequilibrium/ld_decay_40k.hap.ld",h=TRUE)
head(ld)
dim(ld)

# Compute distance between SNP pairs
ld$dist <- abs(ld$POS2 - ld$POS1)
ld$Dprime_abs <- abs(ld$Dprime)

ld_summary_I <- subset(ld,CHR=="I") %>%
  mutate(bin = cut(dist, breaks=seq(0, max(ld$dist), by=1000))) %>%
  group_by(bin) %>%
  summarise(mean_r2 = mean(Dprime_abs, na.rm=TRUE), mid = mean(dist, na.rm=TRUE))

ld_summary_II <- subset(ld,CHR=="II") %>%
  mutate(bin = cut(dist, breaks=seq(0, max(ld$dist), by=1000))) %>%
  group_by(bin) %>%
  summarise(mean_r2 = mean(Dprime_abs, na.rm=TRUE), mid = mean(dist, na.rm=TRUE))

ld_summary_III <- subset(ld,CHR=="III") %>%
  mutate(bin = cut(dist, breaks=seq(0, max(ld$dist), by=1000))) %>%
  group_by(bin) %>%
  summarise(mean_r2 = mean(Dprime_abs, na.rm=TRUE), mid = mean(dist, na.rm=TRUE))

plot(ld_summary_III$mid,ld_summary_III$mean_r2,type='l',ylim=c(0,1),col="red")
lines(ld_summary_II$mid,ld_summary_II$mean_r2,type='l',col="blue")
lines(ld_summary_I$mid,ld_summary_I$mean_r2,type='l',col="magenta")


### R2 ?
ld_summary_I <- subset(ld,CHR=="I") %>%
  mutate(bin = cut(dist, breaks=seq(0, max(ld$dist), by=1000))) %>%
  group_by(bin) %>%
  summarise(mean_r2 = mean(R.2, na.rm=TRUE), mid = mean(dist, na.rm=TRUE))

ld_summary_II <- subset(ld,CHR=="II") %>%
  mutate(bin = cut(dist, breaks=seq(0, max(ld$dist), by=1000))) %>%
  group_by(bin) %>%
  summarise(mean_r2 = mean(R.2, na.rm=TRUE), mid = mean(dist, na.rm=TRUE))

ld_summary_III <- subset(ld,CHR=="III") %>%
  mutate(bin = cut(dist, breaks=seq(0, max(ld$dist), by=1000))) %>%
  group_by(bin) %>%
  summarise(mean_r2 = mean(R.2, na.rm=TRUE), mid = mean(dist, na.rm=TRUE))

ld_summary_IV <- subset(ld,CHR=="IV") %>%
  mutate(bin = cut(dist, breaks=seq(0, max(ld$dist), by=1000))) %>%
  group_by(bin) %>%
  summarise(mean_r2 = mean(R.2, na.rm=TRUE), mid = mean(dist, na.rm=TRUE))

plot(ld_summary_III$mid,ld_summary_III$mean_r2,type='l',ylim=c(0,1),col="red")
lines(ld_summary_II$mid,ld_summary_II$mean_r2,type='l',col="blue")
lines(ld_summary_I$mid,ld_summary_I$mean_r2,type='l',col="magenta")
lines(ld_summary_IV$mid,ld_summary_IV$mean_r2,type='l',col="darkgreen")

ld_summary_I.pos <- subset(ld,CHR=="I" & dist<1000) %>%
  mutate(bin = cut(POS1, breaks=seq(0, max(ld$POS1), by=100000))) %>%
  group_by(bin) %>%
  summarise(mean_r2 = mean(R.2, na.rm=TRUE), mid = mean(POS1, na.rm=TRUE))

ld_summary_I.pos.sup <- subset(ld,CHR=="I" & dist>=1000 & dist<3000) %>%
  mutate(bin = cut(POS1, breaks=seq(0, max(ld$POS1), by=100000))) %>%
  group_by(bin) %>%
  summarise(mean_r2 = mean(R.2, na.rm=TRUE), mid = mean(POS1, na.rm=TRUE))

ld_summary_I.pos.sup10 <- subset(ld,CHR=="I" & dist>=10000) %>%
  mutate(bin = cut(POS1, breaks=seq(0, max(ld$POS1), by=100000))) %>%
  group_by(bin) %>%
  summarise(mean_r2 = mean(R.2, na.rm=TRUE), mid = mean(POS1, na.rm=TRUE))

ld_summary_I.pos.sup30 <- subset(ld,CHR=="I" & dist>=30000) %>%
  mutate(bin = cut(POS1, breaks=seq(0, max(ld$POS1), by=100000))) %>%
  group_by(bin) %>%
  summarise(mean_r2 = mean(R.2, na.rm=TRUE), mid = mean(POS1, na.rm=TRUE))

ld_summary_II.pos <- subset(ld,CHR=="II" & dist<1000) %>%
  mutate(bin = cut(POS1, breaks=seq(0, max(ld$POS1), by=100000))) %>%
  group_by(bin) %>%
  summarise(mean_r2 = mean(R.2, na.rm=TRUE), mid = mean(POS1, na.rm=TRUE))

ld_summary_III.pos <- subset(ld,CHR=="III" & dist<1000) %>%
  mutate(bin = cut(POS1, breaks=seq(0, max(ld$POS1), by=100000))) %>%
  group_by(bin) %>%
  summarise(mean_r2 = mean(R.2, na.rm=TRUE), mid = mean(POS1, na.rm=TRUE))

plot(ld_summary_I.pos$mid,ld_summary_I.pos$mean_r2,type='l',col="red",ylim=c(0,1))
lines(ld_summary_I.pos.sup$mid,ld_summary_I.pos.sup$mean_r2,col="red",lty=2)
lines(ld_summary_I.pos.sup10$mid,ld_summary_I.pos.sup10$mean_r2,col="black",lty=1)
lines(ld_summary_I.pos.sup30$mid,ld_summary_I.pos.sup30$mean_r2,col="blue",lty=1)

plot(ld_summary_II.pos$mid,ld_summary_II.pos$mean_r2,type='l',col="blue")
plot(ld_summary_III.pos$mid,ld_summary_III.pos$mean_r2,type='l',col="magenta")


