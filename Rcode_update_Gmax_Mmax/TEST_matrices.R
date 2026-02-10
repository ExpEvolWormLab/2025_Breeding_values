

Gmat_NaCl <- fit.NaCl$sigma$`u:line`
Rmat_NaCl <- fit.NaCl$sigma$`u:units`

Gmat_NGM <- fit_NGM$sigma$`u:line`
Rmat_NGM <- fit_NGM$sigma$`u:units`

abs(eigen(Gmat_NaCl)$values)/sum(abs(eigen(Gmat_NaCl)$values))
abs(eigen(Gmat_NGM)$values)/sum(abs(eigen(Gmat_NGM)$values))

barplot(eigen(Gmat_NGM)$values)

eigen(Gmat_NGM)$vector
eigen(Gmat_NaCl)$vector

eigen(cov(phen_NaCl_sub[,2:8]))$vector

cor(eigen(Gmat_NGM)$vector,eigen(Gmat_NaCl)$vector)

cor(beta,eigen(Gmat_NGM)$vector)
cor(beta,eigen(Gmat_NaCl)$vector)

cor(eigen(gamma)$vectors[,7],eigen(Gmat_NGM)$vector)
cor(eigen(gamma)$vectors[,7],eigen(Gmat_NaCl)$vector)



gamma