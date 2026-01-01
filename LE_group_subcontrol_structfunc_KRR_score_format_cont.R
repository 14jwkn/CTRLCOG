# For a given k for clustering, type of SC normalization, and type and 
# percentage of thresholding, collect decomposed accuracy scores on the model 
# predicting each cognitive variable in the given batch using controllability 
# calculated with the continuous time system calculated from the set of state 
# matrices specified. Reformat to generate percentages of the total accuracy 
# score.
# Output:
# testacc_noP.csv Test set CV accuracy score with no p-values.
# csplitver,'.csv' Decomposed CV accuracy scores.
# csplitver,'_perc.csv' Decomposed CV accuracy scores with percentages.

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

#Set base path and predict path.
subgroup <- 'full'
sc_subgroup <- 'dr_full'
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/SC_dFC/',
                   sc_subgroup,'/collect/',threstype,'/',thresval,'/')
krrpath <- paste0(basepath,'KRRXFS/cont_',controltype,'_',statetype,'_',septype,'_',
                   nrep,'_',inner_k,'_',outer_k,'_',sctype,'/')

#Set parameters.
nk <- as.numeric(k)
if (septype=='comCFAng') {
  coglist <- c('gCFA','P24_CR','PV','gFngCFA','gCngCFA')
}
ncog <- length(coglist)
if (statetype == 'SC_sFC_dFCcat') {
  states <- c('sc','sFC',paste0('s',1:nk))
} else if (statetype == 'SC_dFCcat') {
  states <- c('sc',paste0('s',1:nk))
} else if (statetype == 'dFCcat') {
  states <- paste0('s',1:nk)
} else if (statetype == 'SC') {
  states <- c('sc')
} else if (statetype == 'sFC') {
  states <- c('sFC')
}
nspace <- length(feature_names)

# Accuracy ----------------------------------------------------------------

#Read in the test accuracy.
infile <- paste0(krrpath,'score_collect.h5')
store <- pd$HDFStore(infile,'r')
inkey <- '/rep_testacc'
testacc <- store$select(inkey)
store$close()
testacc_labs <- colnames(testacc)
testacc <- colMeans(testacc)
testacc <- data.frame(testacc)
colnames(testacc) <- c('Accuracy')

#Save.
outfile <- paste0(krrpath,'testacc_noP.csv')
write.csv(testacc,outfile,row.names=T)

#Read in split accuracy versions.
splitvers <- c('splitacc','splitacc_strict')
nsplitvers <- length(splitvers)
for (spidx in 1:nsplitvers) {
  
  #Read.
  csplitver <- splitvers[spidx]
  store <- pd$HDFStore(infile,'r')
  inkey <- paste0('/rep_',csplitver)
  csplitmat <- store$select(inkey)
  store$close()
  csplitmat <- colMeans(csplitmat)
  splitacc <- array(csplitmat,dim=c(nspace,ncog))
  splitacc <- rbind(splitacc,colSums(splitacc))
  rownames(splitacc) <- c(feature_names,'Total')
  colnames(splitacc) <- coglist
  
  #Save.
  outfile <- paste0(krrpath,csplitver,'.csv')
  write.csv(splitacc,outfile,row.names=T)
  
  #Generate percentage.
  splitacc_perc <- matrix(NA,nrow=nrow(splitacc),ncol=ncol(splitacc))
  rownames(splitacc_perc) <- rownames(splitacc)
  colnames(splitacc_perc) <- colnames(splitacc)
  ctot <- splitacc['Total',]
  for (cfeat in feature_names) {
    crow <- splitacc[cfeat,]
    splitacc_perc[cfeat,] <- crow/ctot
  }
  splitacc_perc <- splitacc_perc*100
  splitacc_perc['Total',] <- splitacc['Total',]
  
  #Save.
  outfile <- paste0(krrpath,csplitver,'_perc.csv')
  write.csv(splitacc_perc,outfile,row.names=T)
}
print('Saved scores.')
