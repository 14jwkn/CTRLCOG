# For a given k for clustering, type of SC normalization, and type and percentage 
# of thresholding, collect feature scores from MST and no MST versions predicting
# each cognitive variable and correlate them.
# Output:
# covha_',ccog,'_pr_mst.csv Correlations between MST and no MST feature scores for each state for one cognitive variable.

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

#Set MST path and normal path.
subgroup <- 'full'
sc_subgroup <- 'dr_full'
mstpath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                  subgroup,'/',k,'/SC_dFC/',sc_subgroup,'/collect/',threstype,'_mst/',
                  thresval,'/KRRXFS/',controltype,'_',statetype,'_',
                  septype,'_',nrep,'_',inner_k,'_',outer_k,'_',sctype)
nopath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                 subgroup,'/',k,'/SC_dFC/',sc_subgroup,'/collect/',threstype,'/',
                 thresval,'/KRRXFS/',controltype,'_',statetype,'_',
                 septype,'_',nrep,'_',inner_k,'_',outer_k,'_',sctype)
  
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

  #Read in the MST scores.
  infile <- paste0(mstpath,'/score_collect.h5')
  inkey <- paste0('/rep_covha_feat_',ccog)
  store <-  pd$HDFStore(infile,'r')
  fold_mst <- store$select(inkey)
  store$close()
  mean_mst <- colMeans(fold_mst)

  #Read in the normal scores.
  infile <- paste0(nopath,'/score_collect.h5')
  inkey <- paste0('/rep_covha_feat_',ccog)
  store <-  pd$HDFStore(infile,'r')
  fold_no <- store$select(inkey)
  store$close()
  mean_no <- colMeans(fold_no)
  
  #Find intersection, feature labels that don't appear in normal are removed.
  featlabs <- colnames(fold_mst)[colnames(fold_mst) %in% colnames(fold_no)]
  mean_mst <- mean_mst[featlabs]
  fold_mst <- fold_mst[,featlabs]
  mean_no <- mean_no[featlabs]
  fold_no <- fold_no[,featlabs]
  
  #Correlate features for each state for fold mean values.
  mean_corr <- matrix(NA,1,(nstates))
  rownames(mean_corr) <- c('mean')
  colnames(mean_corr) <- states
  for (cstate in states) {
    wantlabs <- featlabs[grepl(paste0(cstate,'_'),featlabs)]
    mean_corr['mean',cstate] <- cor(mean_mst[wantlabs],mean_no[wantlabs],method='pearson')
  }
  
  #Save the combined dataframe.
  outmat <- mean_corr
  outfile <- paste0(mstpath,'/covha_',ccog,'_pr_mst.csv')
  write.csv(outmat,outfile,row.names=T)
}
