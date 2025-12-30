


a=read.table("/Users/mallard/Documents/GitHub/2025_Breeding_values/SliM_simulations/slimSimul1G0.txt",h=TRUE,sep='\t')
b=read.table("/Users/mallard/Documents/GitHub/2025_Breeding_values/SliM_simulations/slimSimul1G36.txt",h=TRUE,sep='\t')
c=read.table("/Users/mallard/Documents/GitHub/2025_Breeding_values/SliM_simulations/slimSimul1G50.txt",h=TRUE,sep='\t')

names(a)[2]="freq0"
names(b)[2]="freq36"
names(c)[2]="freq50"

d=merge(a,b)
d=merge(d,c)

head(d)

library("poolSeq")

a0=c();a36=c();a50=c()
for(i in 1:nrow(d)){
  a0=c(a0,  rbinom(1,50,d[i,2]))
  a36=c(a36,  rbinom(1,50,d[i,3]))
  a50=c(a50,  rbinom(1,50,d[i,4]))
}

test36 <- chi.sq.test(a0,50-a0,a36,50-a36)
plot(d[,1],-log10(p.adjust(test36,method='BH')),pch=16,cex=.6)

test50 <- chi.sq.test(a0,50-a0,a50,50-a50)
plot(d[,1],-log10(p.adjust(test36,method='BH')),pch=16,cex=.6)



