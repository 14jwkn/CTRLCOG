# For a given k for clustering, type of SC normalization, and type and percentage 
# of thresholding, find regions with Gramian errors and plot them.
# Output:
# cstate,'_err.jpg Brain plot of regions with Gramian errors.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggseg)
library(ggsegGlasser)
library(scico)

#Catches arguments.
args = commandArgs(trailingOnly=T)
k <- args[1]
sctype <- args[2]
threstype <- args[3]
thresval <- args[4]
cat(k,sctype,threstype,thresval,'\n')

#Set base path.
subgroup <- 'full'
sc_subgroup <- 'dr_full'
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/SC_dFC/',sc_subgroup,'/collect/',
                   threstype,'/',thresval,'/state_images')

#Set parameters.
xwidth <- 3
xheight <- 2.2

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

#Set states.
nk <- as.numeric(k)
states <- c(paste0(sctype,'_SC'),'sFC',paste0('s',as.character(1:nk)))

#For the possible states, see if there exists an error file.
for (cstate in states) {
  cfile <- paste0(basepath,'/',cstate,'_err.csv')
  if (file.exists(cfile)) {
    
    #Read in the file and set 1 for error regions.
    inmat <- read.csv(cfile)
    vmat <- matrix(0,nrow=360,ncol=1)
    rownames(vmat) <- as.character(c(1:360))
    vmat[inmat[,'Region'],] <- 1
    colnames(vmat) <- 'Score'
    
    #Bind atlas and brain image labels.
    curr_v <- cbind(atlasord,vmat)
    atlas_v <- inner_join(base_atlas,curr_v)
    
    #Plot.
    p1 <- atlas_v %>%
      ggplot() +
      geom_brain(atlas=glasser, 
                 position=position_brain(hemi~side),
                 aes(fill=Score)) +
      scale_fill_scico(palette='berlin',
                       midpoint=0) +
      theme_classic() +
      theme(plot.title = element_text(hjust=0.05),
            axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position='none',
            title=element_blank(),
            legend.margin = margin(t=0,
                                   r=17, 
                                   b=0, 
                                   l=0, unit = "pt"),
            plot.margin=grid::unit(c(0,0,0,0),'mm')) 

  
    #Save.
    outfile <- paste0(basepath,'/',cstate,'_err.jpg')
    ggsave(outfile,p1,units='in',width=xwidth,height=xheight,dpi=360)
  }
}
