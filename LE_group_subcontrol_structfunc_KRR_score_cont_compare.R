# For a given k for clustering, type of SC normalization, type and percentage 
# of thresholding, and control metric, collect feature scores from discrete and 
# continuous time system versions predicting each cognitive variable and correlate them.
# Output:
# covha_',ccog,'_pr_divco.csv Correlations between discrete and continuous feature scores for each state for one cognitive variable.

#Set libraries.
library(tidyverse)
library(reticulate)
library(ggplot2)
library(gridExtra)
library(ggseg)
library(ggsegGlasser)
library(scico)
library(ggradar)

#Set virtual environment.
use_virtualenv('../envs/FlexEnv')

#Load python packages.
h5py <- import('h5py')
np <- import('numpy')
pd <- import('pandas')

#Catches arguments.
args = commandArgs(trailingOnly=T)
k <- args[1]
sctype <- args[2]
threstype <- args[3]
thresval <- args[4]
controltype <- args[5]
statetype <- args[6]
septype <- args[7]
nrep <- args[8]
inner_k <- args[9]
outer_k <- args[10]
print(paste(k,sctype,threstype,thresval,controltype,statetype,septype,nrep,
            inner_k,outer_k,sep=' '))

#Set discrete path and continuous path.
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   clean,'/',order,'/',roiname,'/',subgroup,'/',k,'/SC_dFC/',
                   sc_subgroup,'/collect/',threstype,'/',thresval)
discpath <- paste0(basepath,'/KRRXFS/',controltype,'_',statetype,'_',septype,'_',
                   contype,'_',alph,'_',mykernel,'_',ynorm,'_',nrep,'_',inner_k,
                   '_',outer_k,'_',sctype)
contpath <- paste0(basepath,'/KRRXFS/cont_',controltype,'_',statetype,'_',septype,'_',
                   contype,'_',alph,'_',mykernel,'_',ynorm,'_',nrep,'_',inner_k,
                   '_',outer_k,'_',sctype)

#Set parameters.
theme_set(theme_minimal())
nk <- as.numeric(k)
if (septype=='comCFAng') {
  coglist <- c('gCFA','P24_CR','PV','gFngCFA','gCngCFA')
  coglabs <- c('g','gF','gC','gF-g','gC-g')
}
ncog <- length(coglist)
if (statetype == 'SC_dFCcat') {
  states <- c('sc',paste0('s',1:nk))
  statelabs <- c('SC',paste0('S',1:nk))
  feature_names <- c('sc','dFC')
} 
nstates <- length(states)
nspace <- length(feature_names)

# Feature Correlations ----------------------------------------------------

#Cognitive versions.
for (cidx in 1:ncog) {
  ccog <- coglist[cidx]
  ccog_lab <- coglabs[cidx]

  #Read in the discrete scores.
  infile <- paste0(discpath,'/score_collect.h5')
  inkey <- paste0('/rep_',cfeatver,'_feat_',ccog)
  store <-  pd$HDFStore(infile,'r')
  fold_disc <- store$select(inkey)
  store$close()
  mean_disc <- colMeans(fold_disc)

  #Read in the continuous scores.
  infile <- paste0(contpath,'/score_collect.h5')
  inkey <- paste0('/rep_',cfeatver,'_feat_',ccog)
  store <-  pd$HDFStore(infile,'r')
  fold_cont <- store$select(inkey)
  store$close()
  mean_cont <- colMeans(fold_cont)
  
  #Correlate all the features for each state for fold mean values.
  featlabs <- colnames(fold_disc)
  mean_corr <- matrix(NA,1,(nstates))
  rownames(mean_corr) <- c('mean')
  colnames(mean_corr) <- states
  for (cstate in states) {
    wantlabs <- featlabs[grepl(paste0(cstate,'_'),featlabs)]
    mean_corr['mean',cstate] <- cor(mean_disc[wantlabs],mean_cont[wantlabs],method='pearson')
  }
  
  #Save the combined dataframe.
  outmat <- mean_corr
  outfile <- paste0(contpath,'/covha_',ccog,'_pr_divco.csv')
  write.csv(outmat,outfile,row.names=T)
}
