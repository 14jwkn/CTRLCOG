# For a given number of repetitions, inner loop number, and outer loop number for CV, 
# plot the CV average test set accuracy score for different models to illustrate changes
# from varying the threshold type and value, k for clustering, and how the cognitive 
# variable is measured.
# Output:
# 'thresvar_',cctrl,'_',cstate,'_',ck,'.jpg' Plot that illustrates effects of varying the threshold.
# 'kvar_',cctrl,'_',cstate,'_',cthrestype,'_',cthresval,'_',csctype,'.jpg' Plot that illustrates effects of varying k for clustering.
# 'cogvar_',cstate,'_',cthrestype,'_',cthresval,'_',csctype,'_',ck,'.jpg' Plot that illustrates effects of varying the definition of the cognitive variable.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggrepel)

#Catches arguments.
args = commandArgs(trailingOnly=T)
nrep <- args[1]
inner_k <- args[2]
outer_k <- args[3]
cat(nrep,inner_k,outer_k,'\n')

#Set paths.
outpath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                 subgroup,'/multi_param/',sc_subgroup,'/KRRXFS/',
                 nrep,'_',inner_k,'_',outer_k,'/')

#Extract accuracy scores and format to long format, keep values as numbers.
clab <- 'outeracc_allver'
infile <- paste0(outpath,clab,'.csv')
fullmat <- read.csv(infile,row.names=1,header=1)
fullmat <- fullmat[-c((nrow(fullmat)-1):nrow(fullmat)),]
cmat <- fullmat %>%
  rownames_to_column() %>%
  separate(col='rowname',into=c('SC_Type','Control','State','K','Threshold_Type','Threshold_Value'),sep='\\.') %>%
  pivot_longer(cols=-c('SC_Type','Control','State','K','Threshold_Type','Threshold_Value'),
               names_to='Cognition',values_to='Value') %>%
  mutate(Control=case_when(Control=='ave'~'AC',
                           Control=='mod'~'MC')) %>%
  mutate(SC_Type=case_when(SC_Type=='fpt'~'FPT',
                           SC_Type=='raw'~'RAW',
                           SC_Type=='none'~'NONE'))
cmat$Value <- as.numeric(cmat$Value)

#Set plotting path.
plotpath <- paste0(outpath,'specific_plots/')
dir.create(plotpath,recursive=T)

#Set plot parameters.
theme_set(theme_minimal())

# Threshold Variation -----------------------------------------------------

#Set plot parameters.
cwidth <- 10
cheight <- 5
minrange <- -0.1
maxrange <- 0.2

#Select.
ctrltypes <- c('AC','MC')
nctrl <- length(ctrltypes)
cstate <- 'SC_dFCcat'
cstate_lab <- 'SC & dFC'
ck <- '6'
threstypes <- c('none','groupconsist','groupstr')
coglist <- c('g','gCFA','P24_CR','PV','gFngCFA','gCngCFA')
coglabs <- c('g','gF','gC','gF-g','gC-g')

#For each control type.
for (ctrlidx in 1:nctrl) {
  cctrl <- ctrltypes[ctrlidx]
  
  #Extract.
  plotmat <- cmat
  plotmat <- plotmat[plotmat[,'Control']==cctrl,]
  plotmat <- plotmat[plotmat[,'State']==cstate,]
  plotmat <- plotmat[plotmat[,'K']==ck,]
  plotmat <- plotmat[unlist(plotmat[,'Threshold_Type']) %in% threstypes,]
  plotmat <- plotmat[unlist(plotmat[,'Cognition']) %in% coglist,]
  
  #Merge normalization and threshold type and threshold value columns, convert to factors to keep ordering.
  plotmat <- plotmat %>%
    mutate(Category=paste(SC_Type,Threshold_Value,Threshold_Type))
  rawlist <- unique(plotmat$Category)
  sortlist <- c(rawlist[grep('none',rawlist)],
                rawlist[grep('groupconsist',rawlist)],
                rawlist[grep('groupstr',rawlist)])
  plotmat$Category <- factor(plotmat$Category,levels=sortlist)
  plotmat$Cognition <- factor(plotmat$Cognition,levels=coglist)
  
  #Plot.
  p1 <- ggplot(plotmat,aes(x=Category,y=Value,color=Cognition,group=interaction(SC_Type,Threshold_Type,Cognition))) +
    geom_line(linewidth=0.5) + 
    geom_point(size=1) + 
    scale_x_discrete(guide=guide_axis(angle=90)) + 
    scale_y_continuous(breaks=seq(minrange,maxrange,by=0.02)) +
    scale_color_discrete(labels=coglabs) +
    labs(color='Cognition',
         title=cctrl) +
    theme(axis.text.x=element_text(angle=90),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  
  #Save.
  outfile <- paste0(plotpath,'thresvar_',cctrl,'_',cstate,'_',ck,'.jpg')
  ggsave(outfile,p1,units='in',width=cwidth,height=cheight,bg='white')
}

# K Variation -------------------------------------------------------------

#Set plot parameters.
cwidth <- 6
cheight <- 4
minrange <- -0.1
maxrange <- 0.2

#Select.
ctrltypes <- c('AC','MC')
nctrl <- length(ctrltypes)
cstate <- 'SC_dFCcat'
cthrestype <- 'groupconsist'
cthresval <- '50'
csctype <- 'FPT'
coglist <- c('g','gCFA','P24_CR','PV','gFngCFA','gCngCFA')
coglabs <- c('g','gF','gC','gF-g','gC-g')
klist <- as.character(c(2:12))

#For each control type.
for (ctrlidx in 1:nctrl) {
  cctrl <- ctrltypes[ctrlidx]
  
  #Extract.
  plotmat <- cmat
  plotmat <- plotmat[plotmat[,'Control']==cctrl,]
  plotmat <- plotmat[plotmat[,'State']==cstate,]
  plotmat <- plotmat[plotmat[,'Threshold_Type']==cthrestype,]
  plotmat <- plotmat[plotmat[,'Threshold_Value']==cthresval,]
  plotmat <- plotmat[plotmat[,'SC_Type']==csctype,]
  plotmat <- plotmat[unlist(plotmat[,'Cognition']) %in% coglist,]
  
  #Convert to factor to keep ordering.
  plotmat$K <- factor(plotmat$K,levels=klist)
  plotmat$Cognition <- factor(plotmat$Cognition,levels=coglist)
  
  #Plot.
  p1 <- ggplot(plotmat,aes(x=K,y=Value,color=Cognition,group=Cognition)) +
    geom_line(linewidth=0.5) + 
    geom_point(size=1) + 
    scale_x_discrete(guide=guide_axis(angle=90)) + 
    scale_y_continuous(breaks=seq(minrange,maxrange,by=0.02)) +
    scale_color_discrete(labels=coglabs) +
    labs(color='Cognition',
         title=cctrl) +
    theme(axis.text.x=element_text(angle=90),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  
  #Save.
  outfile <- paste0(plotpath,'kvar_',cctrl,'_',cstate,'_',cthrestype,'_',cthresval,'_',csctype,'.jpg')
  ggsave(outfile,p1,units='in',width=cwidth,height=cheight,bg='white')
}

# Cognitive Variation -----------------------------------------------------

#Set plot parameters.
cwidth <- 4
cheight <- 3
minrange <- -0.1
maxrange <- 0.2

#Select.
ctrltypes <- c('AC','MC')
cstate <- 'SC_dFCcat'
ck <- '6'
cthrestype <- 'groupconsist'
cthresval <- '50'
csctype <- 'FPT'
coglist <- c('gCFA','gEFA','gPCA2','PV','NC','P24_CR','NF')
coglabs <- c('g','g EFA','g PCA','gC','gC NIH','gF','gF NIH')
klist <- as.character(c(2:12))
  
#Extract.
plotmat <- cmat
plotmat <- plotmat[unlist(plotmat[,'Control']) %in% ctrltypes,]
plotmat <- plotmat[plotmat[,'State']==cstate,]
plotmat <- plotmat[plotmat[,'K']==ck,]
plotmat <- plotmat[plotmat[,'Threshold_Type']==cthrestype,]
plotmat <- plotmat[plotmat[,'Threshold_Value']==cthresval,]
plotmat <- plotmat[plotmat[,'SC_Type']==csctype,]
plotmat <- plotmat[unlist(plotmat[,'Cognition']) %in% coglist,]

#Convert to factor to keep ordering.
plotmat$K <- factor(plotmat$K,levels=klist)
plotmat$Cognition <- factor(plotmat$Cognition,levels=coglist)
  
#Plot.
p1 <- ggplot(plotmat,aes(x=Cognition,y=Value,color=Control,group=Control)) +
  geom_line(linewidth=0.5) + 
  geom_point(size=1) + 
  scale_x_discrete(guide=guide_axis(angle=90),labels=coglabs) + 
  scale_y_continuous(breaks=seq(minrange,maxrange,by=0.02)) +
  labs(color='Control') +
  theme(axis.text.x=element_text(angle=90),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
  
#Save.
outfile <- paste0(plotpath,'cogvar_',cstate,'_',cthrestype,'_',cthresval,'_',csctype,'_',ck,'.jpg')
ggsave(outfile,p1,units='in',width=cwidth,height=cheight,bg='white')
