## 22_Test_SBxArea_fitness_interaction.R
## Exploratory nonlinear phenotypic selection: does the association of SB
## with relative fitness depend on area.F, conditional on other 5 phenotypes?
## Distinct from locus-by-locus genetic epistasis and from the derived beta.
## Uses exactly the standardized RIL-level fitness data saved in model_beta.

prep_file <- 'RData/Preparatory_GLMM_NaCl_only_withBeta.RData'
out_dir <- 'selection_results/phenotypic_selection_interactions'
traits <- c('SF','SB','FS','FB','BS','BF','area.F')
n_boot <- 2000L
set.seed(20261029)
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(prep_file)&&file.exists(paste0(prep_file,'.RData')))
  prep_file<-paste0(prep_file,'.RData')
e <- new.env(parent=emptyenv());load(prep_file,envir=e)
if(!'model_beta'%in%ls(e))stop('Preparation does not include model_beta; use recomputed script 00')
d <- model.frame(e$model_beta)
if(!all(c('rel_fit',traits)%in%names(d)))
  stop('Saved model_beta does not expose relative fitness and seven traits')
d <- d[complete.cases(d[,c('rel_fit',traits)]),,drop=FALSE]
if(nrow(d)<30L)stop('Insufficient complete fitness-matched RIL observations')

baseline <- lm(reformulate(traits,response='rel_fit'),data=d)
quadratic <- update(baseline, . ~ . + I(SB^2) + I(area.F^2))
interaction <- update(quadratic, . ~ . + SB:area.F)
comparison <- anova(baseline,quadratic,interaction)
write.csv(as.data.frame(comparison),file.path(out_dir,'22_Nested_fitness_models.csv'))

## HC3 heteroskedasticity-robust SE, implemented in base R to avoid dependencies.
HC3 <- function(fit) {
  X<-model.matrix(fit);r<-residuals(fit); h<-hatvalues(fit)
  inv<-solve(crossprod(X))
  v<-inv%*%crossprod(X,X * ((r/pmax(1e-8,1-h))^2))%*%inv
  se<-sqrt(diag(v))
  b<-coef(fit)
  data.frame(term=names(b),estimate=unname(b),SE_HC3=se,
             p_HC3=2*pt(abs(b/se),df.residual(fit),lower.tail=FALSE),
             lower95_HC3=b-qt(.975,df.residual(fit))*se,
             upper95_HC3=b+qt(.975,df.residual(fit))*se,
             row.names=NULL)
}
result<-HC3(interaction)
result$SE_OLS<-coef(summary(interaction))[,2]
result$p_OLS<-coef(summary(interaction))[,4]
write.csv(result,file.path(out_dir,'22_Selection_regression_coefs.csv'),row.names=FALSE)

## Bootstrap whole RIL observations. The original model works on one row
## per RIL; the resample keeps all seven phenotypes and fitness together.
terms <- c('SB','area.F','I(SB^2)','I(area.F^2)','SB:area.F')
boot <- matrix(NA_real_,nrow=n_boot,ncol=length(terms),dimnames=list(NULL,terms))
for(i in seq_len(n_boot)) {
  ix<-sample.int(nrow(d),nrow(d),replace=TRUE)
  f<-try(lm(formula(interaction),data=d[ix,,drop=FALSE]),silent=TRUE)
  if(!inherits(f,'try-error')){
    co<-coef(f)
    boot[i,] <- co[terms]
  }
}
B <- data.frame(term=terms,
  median=vapply(seq_along(terms),function(i)median(boot[,i],na.rm=TRUE),numeric(1)),
  lower95=vapply(seq_along(terms),function(i)quantile(boot[,i],.025,na.rm=TRUE),numeric(1)),
  upper95=vapply(seq_along(terms),function(i)quantile(boot[,i],.975,na.rm=TRUE),numeric(1)),
  finite_draws=colSums(is.finite(boot)))
write.csv(B,file.path(out_dir,'22_RIL_bootstrap_interaction.csv'),row.names=FALSE)

## Conditional phenotypic selection gradient for SB, evaluated at mean SB
## and -1,0,+1 SD of area.F. This is association, not genotype interaction.
vals<-c(-1,0,1)
coef_inter<-coef(interaction)
conditional <- data.frame(area_F_SD=vals,
  d_relative_fitness_d_SB_at_SB0 = unname(coef_inter['SB']+
                                           vals*coef_inter['SB:area.F']))
conditional$lower95_boot<-vapply(vals,function(v)quantile(boot[,'SB']+
                                       v*boot[,'SB:area.F'],.025,na.rm=TRUE),numeric(1))
conditional$upper95_boot<-vapply(vals,function(v)quantile(boot[,'SB']+
                                       v*boot[,'SB:area.F'],.975,na.rm=TRUE),numeric(1))
write.csv(conditional,file.path(out_dir,'22_SB_gradient_by_area.csv'),row.names=FALSE)

pdf(file.path(out_dir,'22_SB_by_area_selection_interaction.pdf'),width=7,height=5)
yrange<-range(c(conditional$lower95_boot,conditional$upper95_boot,0),finite=TRUE)
plot(conditional$area_F_SD,conditional$d_relative_fitness_d_SB_at_SB0,
     type='b',pch=16,lwd=2,ylim=yrange,
     xlab='Body area (standardized)',
     ylab='Conditional SB gradient on relative fitness',bty='n',
     main='Does body area modify selection on SB?')
arrows(conditional$area_F_SD,conditional$lower95_boot,
       conditional$area_F_SD,conditional$upper95_boot,
       code=3,angle=90,length=.07)
abline(h=0,col='grey50',lty=2)
dev.off()
cat('Fitness-matched RILs:',nrow(d),'\n')
cat('SB × area.F estimate:',coef_inter['SB:area.F'],'\n')
cat('Outputs:',out_dir,'\n')
cat('Phenotypic correlational selection is not genetic epistasis.\n')
