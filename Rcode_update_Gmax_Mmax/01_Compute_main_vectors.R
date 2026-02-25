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

BLUE1 <- read.table("RData_BLUEs/BLUES_NGM.txt",sep='\t',h=TRUE)
#BLUE2 <- read.table("RData_BLUEs/BLUES_NaCl.txt",sep='\t',h=TRUE)

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

nb_traits=length(vect_P_traits)
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

trait_means <- colMeans(phen_NGM_with_fitness[, vect_P_traits],na.rm=TRUE)
# Define a grid (adjust resolution if needed)
xseq <- seq(min(phen_NGM_with_fitness$SF), max(phen_NGM_with_fitness$SF), length.out = 50)
yseq <- seq(min(phen_NGM_with_fitness$SB), max(phen_NGM_with_fitness$SB), length.out = 50)

grid <- expand.grid(SF = xseq, SB = yseq)

# Add remaining traits at their mean values
for(i in 3:7){
  grid[vect_P_traits[i]] <- trait_means[i]
}
Z <- as.matrix(grid[, vect_P_traits])
linear_term  <- Z %*% beta
quad_term    <- 0.5 * rowSums((Z %*% gamma) * Z)

grid$pred_fitness <- linear_term + quad_term
zmat <- matrix(grid$pred_fitness, nrow = length(xseq), ncol = length(yseq))

persp(xseq, yseq, zmat,
      theta = 45, phi = 30,
      xlab = "SF",
      ylab = "SB",
      zlab = "Fitness",
      col = "lightblue",
      ticktype = "detailed")

library(rgl)
persp3d(xseq, yseq, zmat,
        col = "lightblue",
        xlab = "SF", ylab = "SB", zlab = "Fitness")

library(ggplot2)

ggplot(grid, aes(SF, SB, z = pred_fitness)) +
  geom_contour_filled() +
  labs(title = "Fitness Landscape",
       fill = "Fitness")


ggplot(grid, aes(SF, SB, z = pred_fitness)) +
  geom_contour_filled(alpha = 0.8) +
  geom_point(data = phen_NGM_with_fitness,
             aes(SF, SB),
             color = phen_NGM_with_fitness$rel_fit,
             size = 2,
             alpha = 0.7,inherit.aes = FALSE) +
  labs(title = "Fitness Landscape with Observed Data",
       fill = "Pred. Fitness")



eig <- eigen(gamma)
# eigenvalues and eigenvectors
values <- eig$values
vectors <- eig$vectors


# First eigenvector (largest eigenvalue)
e1 <- vectors[, 1]
# Last eigenvector (smallest eigenvalue)
e7 <- vectors[, length(values)]

z0 <- colMeans(phen_NGM_with_fitness[, vect_P_traits],na.rm=TRUE)

# Range in eigenvector coordinates
s1 <- seq(-3, 3, length.out = 60)
s7 <- seq(-3, 3, length.out = 60)

grid <- expand.grid(s1 = s1, s7 = s7)

Z <- t(apply(grid, 1, function(pt){
  z0 + pt[1] * e1 + pt[2] * e7
}))
Z <- as.matrix(Z)
colnames(Z) <- vect_P_traits

linear_term <- Z %*% beta
quad_term   <- 0.5 * rowSums((Z %*% gamma) * Z)

grid$fitness <- linear_term + quad_term

zmat <- matrix(grid$fitness,
               nrow = length(s1),
               ncol = length(s7))

persp(s1, s7, zmat,
      theta = 40, phi = 30,
      xlab = "Eigenvector 1",
      ylab = "Eigenvector 7",
      zlab = "Fitness",
      col = "lightblue",
      ticktype = "detailed")

library(ggplot2)
ggplot(grid, aes(s1, s7, z = fitness)) +
  geom_contour_filled() +
  labs(x = "Eigenvector 1",
       y = "Eigenvector 7",
       fill = "Fitness",
       title = "Fitness Landscape along γ Eigenvectors")

