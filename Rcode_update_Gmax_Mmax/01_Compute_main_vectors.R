rm(list=ls())
gc()
#load libraries
#library(coda)
#library(matrixStats)
#library(data.table)
#library(lme4)
#library(MCMCglmm)
#library(parallel)


vect_P_traits= c("SF","SB","FS","FB","BS","BF","area")
# Load phenotypic data

BLUE1 <- read.table("BLUES_NGM.txt",sep='\t',h=TRUE)
BLUE2 <- read.table("BLUES_NaCl.txt",sep='\t',h=TRUE)

## Create a data frame

phen_NGM = data.frame(line=BLUE1$pop_label,SF=NA,SB=NA,FS=NA,FB=NA,BS=NA,BF=NA,area=NA)
for(i in 1:nrow(phen_NGM)){
  phen_NGM[i,2:8]=as.numeric(subset(BLUE1,pop_label==phen_NGM$line[i])$BLUE)
}
phen_NGM=unique(phen_NGM)
dim(phen_NGM)

#Load the fertility table from Noble 2017
#load("fertility.rda")
load("data/fertility.rda")
fertility_NGM <- subset(ecoefs,env=="NGM")
phen_NGM_with_fitness = merge(phen_NGM,fertility_NGM[,1:2])
dim(phen_NGM_with_fitness) # 251 measurements

phen_NGM_with_fitness$rel_fit = exp(phen_NGM_with_fitness$fertility)/mean(exp(phen_NGM_with_fitness$fertility),na.rm=TRUE)

for(i in vect_P_traits){
  phen_NGM_with_fitness[,i]=scale(phen_NGM_with_fitness[,i])
}

model_beta <- lm(rel_fit ~ SF+SB+FS+FB+BS+BF+area, data = phen_NGM_with_fitness)
summary(model_beta)

beta <- coef(model_beta)[vect_P_traits]

compute_gamma <- function(phen_df_with_fit,P_traits=c("SF","SB","FS","FB","BS","BF","area")){
model_gamma <- lm(rel_fit ~ (SF+SB+FS+FB+BS+BF+area)^2+
                    I(SF^2)+I(SB^2)+I(FS^2)+I(FB^2)+I(BS^2)+I(BF^2)+I(area^2), data = phen_df_with_fit)
coef_gamma <- coef(model_gamma)
nb_traits = length(P_traits)
gamma <- matrix(0, nrow = nb_traits, ncol = nb_traits, dimnames = list(P_traits, P_traits))
for (i in P_traits) {
  gamma[i, i] <- 2 * coef_gamma[paste0("I(",i, "^2)")]
}
for (i in 1:(nb_traits-1)) {
  for (j in (i+1):nb_traits) {
    term <- paste0(P_traits[i], ":", P_traits[j])
    gamma[P_traits[i], P_traits[j]] <- coef_gamma[term]
    gamma[P_traits[j], P_traits[i]] <- coef_gamma[term]
  }
}

return(gamma)
}

# True gamma 
gamma=compute_gamma(phen_NGM_with_fitness)
barplot(eigen(gamma)$values)

### We need to bootstrap eigenvalues ###
Y = phen_NGM_with_fitness

set.seed(2385439)
B <- 1000
boot_gamma <- replicate(B,
                       compute_gamma(phen_NGM_with_fitness[sample(nrow(phen_NGM_with_fitness), replace = TRUE), ])
)

boot_eigs <- t(apply(boot_gamma,3,function(x){
  eigen(x)$values
}))
head(boot_eigs)
CI.eigs <- t(apply(boot_eigs, 2, quantile, probs = c(0.025, 0.5, 0.975)))

true_eigs <- eigen(gamma, symmetric = TRUE)$values
colnames(CI.eigs) <- c("lower", "median", "upper")

df.eigen <- data.frame(
  eigen = factor(1:length(true_eigs)),
  true = true_eigs,
  lower = CI.eigs[, "lower"],
  upper = CI.eigs[, "upper"]
)
library(ggplot2)

ggplot(df.eigen, aes(x = eigen, y = true)) +
  geom_point(size = 3, color = "red") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  ylab("Eigenvalue (curvature)") +
  xlab("Eigenvector") +
  ggtitle("True eigenvalues with 95% bootstrap confidence intervals")


## True rotation:
## This is the standard technique used in evolutionary quantitative genetics,
## G-matrix evolution, and curvature analyses (e.g., Jones et al., 2012).


eig_true <- eigen(gamma, symmetric = TRUE)
vec_true <- eig_true$vectors
val_true <- eig_true$values

boot_proj <- matrix(NA, nrow = 1000, ncol = nb_traits)

for (b in 1:1000) {
  gamma_b <- boot_gamma[,,b]
  for (i in 1:7) {
    v <- vec_true[, i]
    boot_proj[b, i] <- t(v) %*% gamma_b %*% v
  }
}

CI.eigs.true.rotation <- apply(boot_proj, 2, quantile, probs = c(0.025, 0.5, 0.975))


df.true.CI <- data.frame(axis = factor(1:7),true = val_true,lower = CI.eigs.true.rotation[1, ],upper = CI.eigs.true.rotation[3, ])

ggplot(df.true.CI, aes(axis, true)) +
  geom_point(size = 3, color = "red") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  ylab("Curvature along TRUE eigenvector") +
  xlab("Phenotypic axis (true eigenvectors)") +
  ggtitle("CIs of curvature along fixed phenotypic axes")

vec_true

