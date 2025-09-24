# For the given k for clustering, type of SC normalization, and type and percentage of 
# thresholding, generate box plots for the correlation between average controllability, 
# modal controllability, and degree across subjects and across regions.
# Output:
# threstype,'_',thresval,'_',sctype,'_suball.jpg' Box plots of the correlations across subjects and across regions.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)

#Set arguments.
args <- commandArgs(trailingOnly=T)
k <- args[1]
sctype <- args[2]
threstype <- args[3]
thresval <- args[4]
cat('Doing:',k,sctype,threstype,thresval,'\n')

#Set parameters.
nk <- as.numeric(k)
nroi <- 360
allstates <- c('sc','sFC',paste0('s',as.character(1:nk)))
allstates_labs <- c('SC','sFC',paste0('S',as.character(1:nk)))
nall <- length(allstates)
corrmethod <- 'pearson'
theme_set(theme_minimal())
nwidth <- 5
nheight <- 5
ntxt <- 15

#Set paths.
subgroup <- 'full'
sc_subgroup <- 'dr_full'
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
             subgroup,'/',k,'/SC_dFC/',sc_subgroup,'/collect/',threstype,'/',thresval,'/')
scpath <- paste0('../outputs/d_SC/',sc_subgroup,'/',threstype,'/',thresval,'/',sctype,'/')
sFCpath <- paste0('../outputs/r_sFC/',sc_subgroup,'/',threstype,'/',thresval,'/')
outpath <- paste0(basepath,'corr_diagnostics/')
dir.create(outpath,recursive=T)

# Read --------------------------------------------------------------------

#Read in the full Gramian error table to find states with errors.
infile <- paste0(basepath,'/state_images/',sctype,'_SC_sFC_dFC_gram.csv')
gramall <- read.csv(infile,header=F)
gramstates <- gramall[(gramall[,2]!=0),1]
ngram <- length(gramstates)

#Get dFC AC, MC, and D.
predlabs <- c()
for (kidx in 1:nk) {
  predlabs <- c(predlabs,paste0('s',kidx,'_r',1:nroi))
}
infile <- paste0(basepath,'ave_tab.csv')
dFCave <- read.csv(infile,row.names=1,header=F)
infile <- paste0(basepath,'mod_tab.csv')
dFCmod <- read.csv(infile,row.names=1,header=F)
infile <- paste0(basepath,'deg_tab.csv')
dFCdeg <- read.csv(infile,row.names=1,header=F)
colnames(dFCave) <- predlabs
colnames(dFCmod) <- predlabs
colnames(dFCdeg) <- predlabs

#Remove dFC errors.
if (threstype != 'substr') {
  for (kidx in 1:nk) {
    errstate <- paste0('s',as.character(kidx))
    if (errstate %in% gramstates) {
      
      #Read regions to remove.
      infile <- paste0(basepath,'/state_images/',errstate,'_err.csv')
      errmat <- read.csv(infile,header=T)
      regrem <- paste0(errstate,'_r',errmat[,'Region'])
      keeplabs <- predlabs[!(predlabs %in% regrem)]
      
      #Remove columns.
      dFCave <- dFCave[,keeplabs]
      dFCmod <- dFCmod[,keeplabs]
      dFCdeg <- dFCdeg[,keeplabs]
      predlabs <- keeplabs
    }
  }
}
  
#Get SC AC, MC, and D.
predlabs <- paste0('sc_r',1:nroi)
infile <- paste0(scpath,'orig_ave_sc.csv')
scave <- read.csv(infile,row.names=1,header=F)
infile <- paste0(scpath,'orig_mod_sc.csv')
scmod <- read.csv(infile,row.names=1,header=F)
infile <- paste0(scpath,'deg_sc.csv')
scdeg <- read.csv(infile,row.names=1,header=F)
colnames(scave) <- predlabs
colnames(scmod) <- predlabs
colnames(scdeg) <- predlabs

#Remove SC errors.
if (threstype != 'substr') {
  errstate <- 'SC'
  if (errstate %in% gramstates) {
    
    #Read regions to remove.
    infile <- paste0(basepath,'/state_images/',paste0(sctype,'_',errstate),'_err.csv')
    errmat <- read.csv(infile,header=T)
    regrem <- paste0('sc_r',errmat[,'Region'])
    keeplabs <- predlabs[!(predlabs %in% regrem)]
    
    #Remove columns.
    scave <- scave[,keeplabs]
    scmod <- scmod[,keeplabs]
    scdeg <- scdeg[,keeplabs]
    predlabs <- keeplabs
  }
}

#Get sFC AC, MC, and D.
predlabs <- paste0('sFC_r',1:nroi)
infile <- paste0(sFCpath,'ave_sFC.csv')
sFCave <- read.csv(infile,row.names=1,header=F)
infile <- paste0(sFCpath,'mod_sFC.csv')
sFCmod <- read.csv(infile,row.names=1,header=F)
infile <- paste0(sFCpath,'deg_sFC.csv')
sFCdeg <- read.csv(infile,row.names=1,header=F)
colnames(sFCave) <- predlabs
colnames(sFCmod) <- predlabs
colnames(sFCdeg) <- predlabs

#Remove sFC errors.
if (threstype != 'substr') {
  errstate <- 'sFC'
  if (errstate %in% gramstates) {
    
    #Read regions to remove.
    infile <- paste0(basepath,'/state_images/',errstate,'_err.csv')
    errmat <- read.csv(infile,header=T)
    regrem <- paste0('sFC_r',errmat[,'Region'])
    keeplabs <- predlabs[!(predlabs %in% regrem)]
    
    #Remove columns.
    sFCave <- sFCave[,keeplabs]
    sFCmod <- sFCmod[,keeplabs]
    sFCdeg <- sFCdeg[,keeplabs]
    predlabs <- keeplabs
  }
}

#Concatenate.
allave <- cbind(scave,sFCave,dFCave)
allmod <- cbind(scmod,sFCmod,dFCmod)
alldeg <- cbind(scdeg,sFCdeg,dFCdeg)
nsub <- nrow(allave)
sublist <- rownames(allave)
full_featlabs <- colnames(allave)

# Collect -----------------------------------------------------------------

#Initialize.
tmult <- 1.2
ymult <- 1
cpad <- 10
plist <- list()
pltidx <- 1

# Region-wise -------------------------------------------------------------

#Find AC-D correlations between subjects for every region for each state.
acd_mat <- matrix(NA,nrow=nroi,ncol=nall)
rownames(acd_mat) <- paste0('r',as.character(1:nroi))
colnames(acd_mat) <- allstates
for (cstate in allstates) {
  
  #Extract.
  cstate_labs <- grepl(cstate,full_featlabs,fixed=T)
  cave <- allave[,cstate_labs]
  cdeg <- alldeg[,cstate_labs]
  
  #Strip away the state label.
  creg_labs <- colnames(cdeg)
  creg_labs <- str_remove(creg_labs,paste0(cstate,'_'))
  c_nroi <- length(creg_labs)
  
  #For each region.
  for (ridx in 1:c_nroi) {
    cregion <- creg_labs[ridx]
    cregion_lab <- paste0(cstate,'_',cregion)
    cave_reg <- cave[,cregion_lab]
    cdeg_reg <- cdeg[,cregion_lab]
    acd_mat[cregion,cstate] <- cor(cave_reg,cdeg_reg)
  }
}

#Reformat names.
colnames(acd_mat) <- allstates_labs

#Convert to long format and drop NA and convert to plotting labels.
plotmat <- data.frame(acd_mat) %>%
  pivot_longer(everything()) %>%
  drop_na() %>%
  mutate(name=replace(name,name=='sc','SC'))
plotmat$name <- factor(plotmat$name,levels=allstates_labs)

#Plot boxplots.
p2 <- ggplot(plotmat,aes(x=name,y=value)) +
    geom_boxplot() +
    scale_y_continuous(limits=c(-1,1),breaks=pretty(seq(-1,1,by=0.1),n=6)) +
    ggtitle('AC vs S') +
    ylab('Participant Correlation Per Region') +
    theme(plot.title=element_text(size=(ntxt*tmult),hjust=0.5),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=(ntxt*ymult),margin=margin(t=0,r=cpad,b=0,l=0)),
          panel.grid.minor.y=element_blank(),
          panel.grid.minor.x=element_blank(),
          text=element_text(size=ntxt)) 
plist[[pltidx]] <- p2
pltidx <- pltidx + 1

#Find MC-D correlations between subjects for every region for each state.
mcd_mat <- matrix(NA,nrow=nroi,ncol=nall)
rownames(mcd_mat) <- paste0('r',as.character(1:nroi))
colnames(mcd_mat) <- allstates
for (cstate in allstates) {
  
  #Extract.
  cstate_labs <- grepl(cstate,full_featlabs,fixed=T)
  cmod <- allmod[,cstate_labs]
  cdeg <- alldeg[,cstate_labs]
  
  #Strip away the state label.
  creg_labs <- colnames(cdeg)
  creg_labs <- str_remove(creg_labs,paste0(cstate,'_'))
  c_nroi <- length(creg_labs)
  
  #For each region.
  for (ridx in 1:c_nroi) {
    cregion <- creg_labs[ridx]
    cregion_lab <- paste0(cstate,'_',cregion)
    cmod_reg <- cmod[,cregion_lab]
    cdeg_reg <- cdeg[,cregion_lab]
    mcd_mat[cregion,cstate] <- cor(cmod_reg,cdeg_reg)
  }
}

#Reformat names.
colnames(mcd_mat) <- allstates_labs

#Convert to long format and drop NA and convert to plotting labels.
plotmat <- data.frame(mcd_mat) %>%
  pivot_longer(everything()) %>%
  drop_na() %>%
  mutate(name=replace(name,name=='sc','SC'))
plotmat$name <- factor(plotmat$name,levels=allstates_labs)

#Plot boxplots.
p2 <- ggplot(plotmat,aes(x=name,y=value)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-1,1),breaks=pretty(seq(-1,1,by=0.1),n=6)) +
  ggtitle('MC vs S') +
  ylab('') +
  theme(plot.title=element_text(size=(ntxt*tmult),hjust=0.5),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=(ntxt*ymult),margin=margin(t=0,r=cpad,b=0,l=0)),
        panel.grid.minor.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        text=element_text(size=ntxt)) 
plist[[pltidx]] <- p2
pltidx <- pltidx + 1

#Find AC-MC correlations between subjects for every region for each state.
acmc_mat <- matrix(NA,nrow=nroi,ncol=nall)
rownames(acmc_mat) <- paste0('r',as.character(1:nroi))
colnames(acmc_mat) <- allstates
for (cstate in allstates) {
  
  #Extract.
  cstate_labs <- grepl(cstate,full_featlabs,fixed=T)
  cave <- allave[,cstate_labs]
  cmod <- allmod[,cstate_labs]
  
  #Strip away the state label.
  creg_labs <- colnames(cave)
  creg_labs <- str_remove(creg_labs,paste0(cstate,'_'))
  c_nroi <- length(creg_labs)
  
  #For each region.
  for (ridx in 1:c_nroi) {
    cregion <- creg_labs[ridx]
    cregion_lab <- paste0(cstate,'_',cregion)
    cave_reg <- cave[,cregion_lab]
    cmod_reg <- cmod[,cregion_lab]
    acmc_mat[cregion,cstate] <- cor(cave_reg,cmod_reg)
  }
}

#Reformat names.
colnames(acmc_mat) <- allstates_labs

#Convert to long format and drop NA and convert to plotting labels.
plotmat <- data.frame(acmc_mat) %>%
  pivot_longer(everything()) %>%
  drop_na() 
plotmat$name <- factor(plotmat$name,levels=allstates_labs)

#Plot boxplots.
p2 <- ggplot(plotmat,aes(x=name,y=value)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-1,1),breaks=pretty(seq(-1,1,by=0.1),n=6)) +
  ggtitle('AC vs MC') +
  ylab('') +
  theme(plot.title=element_text(size=(ntxt*tmult),hjust=0.5),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=(ntxt*ymult),margin=margin(t=0,r=cpad,b=0,l=0)),
        panel.grid.minor.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        text=element_text(size=ntxt)) 
plist[[pltidx]] <- p2
pltidx <- pltidx + 1

# Subject-wise ------------------------------------------------------------

#Find AC-D correlations between regions for every subject for each state.
acd_mat <- matrix(NA,nrow=nsub,ncol=nall)
rownames(acd_mat) <- sublist
colnames(acd_mat)  <- allstates
for (cstate in allstates) {
  
  #Extract.
  cstate_labs <- grepl(cstate,full_featlabs,fixed=T)
  cave <- allave[,cstate_labs]
  cdeg <- alldeg[,cstate_labs]
  
  #For every subject.
  for (csub in sublist) {
    cave_sub <- unlist(cave[csub,])
    cdeg_sub <- unlist(cdeg[csub,])
    acd_mat[csub,cstate] <- cor(cave_sub,cdeg_sub)
  }
}

#Reformat names.
colnames(acd_mat) <- allstates_labs

#Convert to long format and drop NA and convert to plotting labels.
plotmat <- data.frame(acd_mat) %>%
  pivot_longer(everything()) %>%
  drop_na() %>%
  mutate(name=replace(name,name=='sc','SC'))
plotmat$name <- factor(plotmat$name,levels=allstates_labs)

#Plot boxplots.
p2 <- ggplot(plotmat,aes(x=name,y=value)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-1,1),breaks=pretty(seq(-1,1,by=0.1),n=6)) +
  ggtitle('') +
  ylab('Region Correlation Per Participant') +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=(ntxt*ymult),margin=margin(t=0,r=cpad,b=0,l=0)),
        panel.grid.minor.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        text=element_text(size=ntxt)) 
plist[[pltidx]] <- p2
pltidx <- pltidx + 1

#Find MC-D correlations between regions for every subject for each state.
mcd_mat <- matrix(NA,nrow=nsub,ncol=nall)
rownames(mcd_mat) <- sublist
colnames(mcd_mat) <- allstates
for (cstate in allstates) {
  
  #Extract.
  cstate_labs <- grepl(cstate,full_featlabs,fixed=T)
  cmod <- allmod[,cstate_labs]
  cdeg <- alldeg[,cstate_labs]
  
  #For every subject.
  for (csub in sublist) {
    cmod_sub <- unlist(cmod[csub,])
    cdeg_sub <- unlist(cdeg[csub,])
    mcd_mat[csub,cstate] <- cor(cmod_sub,cdeg_sub)
  }
}

#Reformat names.
colnames(mcd_mat) <- allstates_labs

#Convert to long format and drop NA and convert to plotting labels.
plotmat <- data.frame(mcd_mat) %>%
  pivot_longer(everything()) %>%
  drop_na() %>%
  mutate(name=replace(name,name=='sc','SC'))
plotmat$name <- factor(plotmat$name,levels=allstates_labs)

#Plot boxplots.
p2 <- ggplot(plotmat,aes(x=name,y=value)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-1,1),breaks=pretty(seq(-1,1,by=0.1),n=6)) +
  ggtitle('') +
  ylab('') +
  theme(plot.title=element_text(size=(ntxt*tmult),hjust=0.5),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=(ntxt*ymult),margin=margin(t=0,r=cpad,b=0,l=0)),
        panel.grid.minor.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        text=element_text(size=ntxt)) 
plist[[pltidx]] <- p2
pltidx <- pltidx + 1

#Find AC-MC correlations between regions for every subject for each state.
acmc_mat <- matrix(NA,nrow=nsub,ncol=nall)
rownames(acmc_mat) <- sublist
colnames(acmc_mat) <- allstates
for (cstate in allstates) {
  
  #Extract.
  cstate_labs <- grepl(cstate,full_featlabs,fixed=T)
  cave <- allave[,cstate_labs]
  cmod <- allmod[,cstate_labs]
  
  #For every subject.
  for (csub in sublist) {
    cave_sub <- unlist(cave[csub,])
    cmod_sub <- unlist(cmod[csub,])
    acmc_mat[csub,cstate] <- cor(cave_sub,cmod_sub)
  }
}

#Reformat names.
colnames(acmc_mat) <- allstates_labs

#Convert to long format and drop NA and convert to plotting labels.
plotmat <- data.frame(acmc_mat) %>%
  pivot_longer(everything()) %>%
  drop_na() %>%
  mutate(name=replace(name,name=='sc','SC'))
plotmat$name <- factor(plotmat$name,levels=allstates_labs)

#Plot boxplots.
p2 <- ggplot(plotmat,aes(x=name,y=value)) +
  geom_boxplot() +
  scale_y_continuous(limits=c(-1,1),breaks=pretty(seq(-1,1,by=0.1),n=6)) +
  ggtitle('') +
  ylab('') +
  theme(plot.title=element_text(size=(ntxt*tmult),hjust=0.5),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=(ntxt*ymult),margin=margin(t=0,r=cpad,b=0,l=0)),
        panel.grid.minor.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        text=element_text(size=ntxt)) 
plist[[pltidx]] <- p2
pltidx <- pltidx + 1

# Collect -----------------------------------------------------------------

#Save.
gatherpath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                    subgroup,'/',k,'/SC_dFC/',
                    sc_subgroup,'/collect/corr_diagnostics/')
dir.create(gatherpath,recursive=T)
outfile <- paste0(gatherpath,threstype,'_',thresval,'_',sctype,'_suball.jpg')
ggsave(outfile,arrangeGrob(grobs=plist,ncol=3),
       units='in',width=(nwidth*3),height=(nheight*2),bg='white')
