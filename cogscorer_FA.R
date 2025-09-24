# For the expanded sample, conduct EFA and CFA and retrieve scores and loadings.
# Adapts code from:
# https://github.com/adolphslab/HCP_MRI-behavior/blob/master/intelligence.ipynb
# Output:
# dubois_efa_load.csv EFA loadings.
# dubous_efa_scores.csv EFA scores.
# dubois_cfa_load.csv CFA loadings.
# dubous_cfa_scores.csv CFA scores.

#Set libraries.
library(tidyverse)
library(psych)
library(lavaan)

#Get cognitive and other predictor data.
subgroup <- 'ourfa_subs'
infile <- paste0('../outputs/c_cognition/',subgroup,'/cogscores.csv')
othercog <- read_csv(infile,col_names=T) %>%
  mutate(Subject=as.character(Subject)) %>%
  mutate_if(is.character,as.factor)

#Isolate desired cognitive variables.
ourcog <- c('PMAT24','PicVoc','ProcSpeed','CardSort','Flanker','ListSort',
            'PicSeq','ReadVoc','WordMem','LineOrient')
cogmat <- as.data.frame(othercog[,ourcog])
rownames(cogmat) <- as.vector(unlist(othercog[,'Subject'],use.names=F))

#Set up helper functions.
# compute Comparative Fit Index for a factor analysis 
CFI <-function(x){
  return((1-((x$STATISTIC-x$dof))/(x$null.chisq-x$null.dof)))
}
# compute Comparative Fit Index for a bifactor analysis 
CFI_biv <-function(x){
  return((1-((x$stats$STATISTIC-x$stats$dof))/(x$stats$null.chisq-x$stats$null.dof)))
}
# compute implied matrix for a factor analysis
impliedMatrix<-function(x){
  if (dim(x$loadings)[2]==1) {
    imp      <- x$loadings %*% t(x$loadings) 
  } else {
    imp      <- x$loadings %*% x$Phi %*% t(x$loadings) 
  }
  diag(imp)<- diag(imp) + x$uniquenesses
  return(imp)
}
# compute implied matrix for a bifactor analysis
impliedMatrix_biv<-function(x){
  Gloadings     <- x$schmid$sl[,1]
  Floadings     <- x$schmid$sl[,2:(ncol(x$schmid$sl)-3)]
  uniquenesses  <- x$schmid$sl[,ncol(x$schmid$sl)-1]
  imp           <- Gloadings %*% t(Gloadings) + Floadings %*% t(Floadings)
  diag(imp)     <- diag(imp) + uniquenesses
  return(imp)
}

#Parallel analysis.
out = fa.parallel(cogdf,plot=F) 
faValues = out$fa.values
faSim = out$fa.sim
faSimR = out$fa.simr

# Dubois EFA --------------------------------------------------------------

#Run exploratory FA.
fm     <- "mle"       # use maximum likelihood estimator
rotate <- "oblimin"   # use oblimin factor rotation
fitInds <- matrix(,nrow=1,ncol=9)
rownames(fitInds) <- c('b4')
colnames(fitInds) <- c('CFI','RMSEA','SRMR','BIC','om_h','om_s1','om_s2','om_s3','om_s4')

#Observed covariance matrices
obs <- cov(cogdf)
lobs <- obs[!lower.tri(obs)]

#Bi-factor model.
model <- 'b4'
b4 <- omega(cogdf,nfactors=4,fm=fm,key=NULL,flip=FALSE,
            digits=3,title="Omega",sl=TRUE,labels=NULL, plot=FALSE,
            n.obs=NA,rotate=rotate,Phi = NULL,option="equal",covar=FALSE)
imp <- impliedMatrix_biv(b4)
limp <- imp[!lower.tri(imp)]
fitInds[model,1] <- CFI_biv(b4)
fitInds[model,2] <- b4$schmid$RMSEA[1]
fitInds[model,3] <- sqrt(mean((limp - lobs)^2))
fitInds[model,4] <- b4$stats$BIC
fitInds[model,5] <- b4$omega_h
fitInds[model,6:9] <- b4$omega.group[-1,3]

#Show output.
print(fitInds,digits=3)
print(b4)
diagram(b4,digits=3,cut=.2)

#Produce loadings.
efa_load <- b4$schmid$sl[,1:5]
colnames(efa_load) <- c('gEFA','cryEFA','spdEFA','visEFA','memEFA')
outfile <- paste0('../outputs/c_cognition/',subgroup,'/dubois_efa_load.csv')
write.csv(efa_load,outfile,row.names=T)

#Produce scores.
efa1_scores <- f1$scores
efa_scores <- factor.scores(cogdf,b4$schmid$sl[,1:5])$scores
colnames(efa_scores) <- c('gEFA','cryEFA','spdEFA','visEFA','memEFA')
efa_scores <- rownames_to_column(data.frame(efa_scores),'Subject')
outfile <- paste0('../outputs/c_cognition/',subgroup,'/dubois_efa_scores.csv')
write.csv(efa_scores,outfile,row.names=F)

# CFA ---------------------------------------------------------------------

#biB enforces loadings of 1 for factors defined by only two observed variables. Model
#does not converge when ListSort is included so did not include it here.
biB <- '
    #g-factor
    g   =~ CardSort + Flanker + ProcSpeed + PicVoc + ReadVoc + PMAT24 + LineOrient + WordMem + PicSeq
    #Domain factors
    spd =~ CardSort + Flanker + ProcSpeed
    cry =~ 1*PicVoc + 1*ReadVoc
    vis =~ 1*PMAT24 + 1*LineOrient   
    mem =~ 1*WordMem + 1*PicSeq
    #Domain factors are not correlated with g
    g ~~ 0*spd
    g ~~ 0*cry
    g ~~ 0*vis
    g ~~ 0*mem
    #Domain factors are not correlated with one another
    spd ~~ 0*cry
    spd ~~ 0*vis
    spd ~~ 0*mem
    cry ~~ 0*vis
    cry ~~ 0*mem
    vis ~~ 0*mem
'
mod_biB <- cfa(biB,data=cogdf,estimator='ML')

#Show output.
print(mod_biB)
print(fitMeasures(mod_biB,c("cfi","tli","rmsea","srmr","aic","bic","chisq","df")))

#Produce loadings.
cfa_load <- data.frame(t(inspect(mod_biB,what='std')$lambda))
cfa_load$ListSort <- rep(0,times=nrow(cfa_load)) #Add back ListSort.
cfa_load <- t(cfa_load[,ourcog])
colnames(cfa_load) <- paste0(colnames(cfa_load),'CFA')
outfile <- paste0('../outputs/c_cognition/',subgroup,'/dubois_cfa_load.csv')
write.csv(cfa_load,outfile,row.names=T)

#Produce scores.
cfa_scores <- lavPredict(mod_biB)
rownames(cfa_scores) <- as.vector(unlist(othercog[,'Subject'],use.names=F))
colnames(cfa_scores) <- paste0(colnames(cfa_scores),'CFA')

#Save the scores.
cfa_scores <- rownames_to_column(data.frame(cfa_scores),'Subject')
outfile <- paste0('../outputs/c_cognition/',subgroup,'/dubois_cfa_scores.csv')
write.csv(cfa_scores,outfile,row.names=F)
