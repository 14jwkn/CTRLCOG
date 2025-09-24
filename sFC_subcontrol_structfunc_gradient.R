# For the given type and percentage of thresholding, plot each gradient for sFC
# on the brain. Generate a file to adjust the sign of the gradient as desired.
# Output:
# sFC_gradients_flip.csv File to change the gradient sign.
# 'sFC_gradients_',as.character(gidx),'.jpg' Single gradient plotted on the brain.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(ggseg)
library(ggsegGlasser)
library(RColorBrewer)

#Catches arguments.
args <- commandArgs(trailingOnly=T)
threstype <- args[1]
thresval <- args[2]
cat(threstype,thresval,'\n')

#Set paths.
basepath <- paste0('../outputs/r_sFC/',clean,'/',order,'/',roiname,'/',subgroup,'/',
                  threstype,'/',thresval,'/')

#Extract region names.
atlasfile <- '../inputs/data/atlas/atlasorder.csv'
atlasord <- read_csv(atlasfile)[,'labelname'] %>%
  mutate_all(sub,pattern='_ROI',replacement='') %>%
  separate(col='labelname',
           sep='_',
           into=c('hemi','region')) %>%
  mutate(hemi=case_when(
    hemi=='L'~'left',
    hemi=='R'~'right')) %>%
  select('region','hemi')
regname <- paste0(atlasord$region,'_',atlasord$hemi)

#Make base atlas.
base_atlas <- as.data.frame(na.omit(cbind(glasser$data$region,glasser$data$hemi)))
colnames(base_atlas) <- c('region', 'hemi')

#Set colors.
gcolors <- list(c('purple','green'),c('blue','red'))

#Read in the first three gradients.
infile <- paste0(basepath,'sFC_gradients.csv')
ingrad <- read.csv(infile,header=F)[,1:3]
ngrad <- 3

#Generate flippers.
outfile <- paste0(basepath,'sFC_gradients_flip.csv')
if (file.exists(outfile)) {
  yflip <- as.vector(t(read.csv(outfile,header=F)))
} else {
  
  #Set up flippers.
  yflip <- rep(F,each=ngrad)
  
  #Generate file.
  outflip <- t(yflip) 
  outflip[outflip == TRUE] = 'T'
  outflip[outflip == FALSE] = 'F'
  write.table(outflip,outfile,sep=',',row.names=F,col.names=F)
}

#For each gradient, plot.
for (gidx in 1:ngrad) {
  cgrad <- ingrad[,gidx]
  ccolors <- gcolors[[gidx]]
  
  #Flip.
  if (yflip[gidx]) {
    cgrad <- -cgrad
  }
  
  #Plot.
  Score <- cgrad
  curr_v <- cbind(atlasord,Score)
  atlas_v <- inner_join(base_atlas,curr_v)
  ggobj <- atlas_v %>%
    ggplot() +
    geom_brain(atlas=glasser, 
               position=position_brain(hemi~side),
               aes(fill=Score)) +
    scale_fill_gradient(low=ccolors[1],high=ccolors[2],na.value='grey') +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.05),
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.margin = margin(t=0,
                                 r=17, 
                                 b=0, 
                                 l=0, unit = "pt"),
          plot.margin=grid::unit(c(1,0,1,0),'mm')) 
  outfile <- paste0(basepath,'sFC_gradients_',as.character(gidx),'.jpg')
  ggsave(outfile,ggobj,units='in',width=(3*2),height=(2.2*2),dpi=720)
}
cat('Saved.\n')
