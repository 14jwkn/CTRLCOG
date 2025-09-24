# For the given k for clustering, type of SC normalization, and type and percentage of 
# thresholding, plot average controllability, modal controllability, and degree
# on the brain along with radar plots to examine regions with the highest values.
# Output:
# threectrl.jpg Unthresholded colored brain plots for average controllability, modal controllability, and degree for SC, sFC, and dFC.
# threectrl_legend.jpg Legends for the unthresholded brain plots.
# cctrl,'_',as.character(cq),'.jpg' Brain plots with thresholding along with radar plots for networks and regions to examine what regions are highest in controllability and degree.

#Set libraries.
library(tidyverse)
library(reticulate)
library(ggplot2)
library(gridExtra)
library(ggseg)
library(ggsegGlasser)
library(scico)
library(RColorBrewer)
library(ggradar)

#Set virtual environment.
use_virtualenv('../envs/FlexEnv')

#Load python packages.
h5py <- import('h5py')
np <- import('numpy')
pd <- import('pandas')

# Functions ---------------------------------------------------------------

#Create label separator, 1 is within label.
labsep <- function(arealabels,arealist) {
  narea <- length(arealist)
  areamat <- matrix(0,nrow=360,ncol=narea)
  colnames(areamat) <- arealist
  for (aidx in 1:narea) {
    areaid <- as.character(aidx)
    areamat[arealabels == areaid,aidx] = 1
  }
  areasum <- t(data.frame(colSums(areamat)))
  return(list(areamat,areasum))
}

#Create dice calculator.
dice <- function(ab,a,b) {
  dice <- 2*(ab/(a+b))
  return (dice)
}

#Create area plotter with radar maximization for only one sign and with no averages (issue with range).
uni_areaplot <- function(plotmat,cstate_lab,arealabels,arealist,nareas,areamat,areasum,areacolors,
                         rad_txtsize,atlasord,base_atlas,nacolor,radcolor,plist,pltidx) {
  
  #Extract loadings for areas.
  cstate_plot <- plotmat
  cstate_plot[plotmat>0] <- 1
  cstate_plot <- cstate_plot*as.numeric(arealabels)
  cstate_plot <- arealist[cstate_plot]
  Score <- factor(cstate_plot,levels=arealist)
  
  #Plot.
  curr_v <- cbind(atlasord,Score)
  atlas_v <- inner_join(base_atlas,curr_v)
  ggobj <- atlas_v %>%
    ggplot() +
    geom_brain(atlas=glasser, 
               position=position_brain(hemi~side),
               aes(fill=Score)) +
    scale_fill_manual(values=areacolors,na.value=nacolor) +
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
          plot.margin=grid::unit(c(1,0,1,0),'mm')) 
  plist[[pltidx]] <- ggobj
  pltidx <- pltidx + 1
  
  #Project to networks.
  pos.state <- plotmat
  pos.state.rep <- matrix(pos.state,nrow=360,ncol=nareas)
  pos.state.values <- pos.state.rep
  pos.state.values[!areamat] <- NA
  
  #Find dice using count of regions of interest, count of regions in the network,
  #and count of the overlap.
  posA <- sum(!is.na(pos.state))
  posB <- areasum
  posAB <-  colSums(!is.na(pos.state.values))
  pos.state.dice <- dice(posAB,posA,posB)
  
  #Maximize and plot.
  statepos <- sum((plotmat>0),na.rm=T)
  ctitle <- paste0(toupper(cstate_lab),' | P=',as.character(statepos))
  pos.state.sep <- pos.state.dice
  cposmax <- max(pos.state.sep)
  badder = data.frame('Group'=c('P'))
  both.sep <- rbind(pos.state.sep)
  colnames(both.sep) <- colnames(areamat)
  both.sep <- cbind(badder,both.sep)
  maxrad <- max(both.sep[,-1])
  if (maxrad == 0) {
    maxrad <- 1
  }
  midrad <- maxrad/2
  minrad <- 0
  ggobj <- ggradar(both.sep,
                   values.radar = c('','',''),
                   grid.min=minrad,grid.mid=midrad,grid.max=maxrad,
                   axis.label.size=rad_txtsize,
                   group.line.width=1,
                   group.point.size=0,
                   group.colours=radcolor,
                   legend.position='none') +
    theme(plot.margin=grid::unit(c(1,0,1,0),'mm'),
          plot.title=element_text(size=titlesize,hjust=0,face='bold')) +
    labs(title=ctitle)
  plist[[pltidx]] <- ggobj
  pltidx <- pltidx + 1
  
  #Return items.
  return(list(plist,pltidx))
}

# Driver ------------------------------------------------------------------

#Catches arguments.
args = commandArgs(trailingOnly=T)
k <- args[1]
sctype <- args[2]
threstype <- args[3]
thresval <- args[4]
cat(k,sctype,threstype,thresval,'\n')

#Set base path and predict path.
subgroup <- 'full'
sc_subgroup <- 'dr_full'
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                   subgroup,'/',k,'/SC_dFC/',sc_subgroup,'/collect/',threstype,'/',thresval,'/')
scpath <- paste0('../outputs/d_SC/',sc_subgroup,'/',threstype,'/',thresval,'/',sctype,'/')
sFCpath <- paste0('../outputs/r_sFC/',sc_subgroup,'/',threstype,'/',thresval,'/')
outpath <- paste0(basepath,'/state_images/')

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

#Create color palette.
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal,qual_col_pals$maxcolors,rownames(qual_col_pals)))

#Divide networks and give colors.
netlabels <- scan('colenetlabels.txt',character(),quote='')
netlist <- c('VIS1','VIS2','SMN',
             'CON','DAN','LAN',
             'FPN','AUD','DMN',
             'PMM','VMM',
             'ORA')
nnet <- length(netlist)
outmat <- labsep(netlabels,netlist)
netmat <- outmat[[1]]
netsum <- outmat[[2]]
netcolors <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
               '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#AEFD6C', '#fff200')
names(netcolors) <- netlist

#Divide regions and give colors.
reglabels <- read.csv('../inputs/data/atlas/area_glasser_atlas.csv')$Area
reglabels <-  case_when(
  reglabels == 'Primary Visual Cortex' ~ 'V1',
  reglabels == 'MT+ Complex and Neighboring Visual Areas' ~ 'MT+',
  reglabels == 'Dorsal Stream Visual Cortex' ~ 'DRStrm',
  reglabels == 'Early Visual Cortex' ~ 'EVis',
  reglabels == 'Ventral Stream Visual Cortex' ~ 'VSStrm',
  reglabels == 'Somatosensory and Motor Cortex' ~ 'SM',
  reglabels == 'Premotor Cortex' ~ 'PreM',
  reglabels == 'Posterior Cingulate Cortex' ~ 'PC',
  reglabels == 'Early Auditory Cortex' ~ 'EAud',
  reglabels == 'Temporo-Parieto-Occipital Junction' ~ 'TPOJ',
  reglabels == 'DorsoLateral Prefrontal Cortex' ~ 'DLPFC',
  reglabels == 'Superior Parietal Cortex' ~ 'SuPar',
  reglabels == 'Paracentral Lobular and Mid Cingulate Cortex'~'ParaCin',
  reglabels == 'Anterior Cingulate and Medial Prefrontal Cortex'~'ACMP',
  reglabels == 'Orbital and Polar Frontal Cortex'~'OPF',
  reglabels == 'Inferior Frontal Cortex'~'IFron',
  reglabels == 'Posterior Opercular Cortex'~'POp',
  reglabels == 'Insular and Frontal Opercular Cortex'~'IFOp',
  reglabels == 'Auditory Association Cortex'~'AA',
  reglabels == 'Inferior Parietal Cortex'~'InPar',
  reglabels == 'Medial Temporal Cortex'~'MTem',
  reglabels == 'Lateral Temporal Cortex'~'LTem')
reglist <- unique(reglabels)
nreg <- length(reglist)
for (ridx in 1:nreg) {
  creg <- reglist[ridx]
  reglabels[reglabels == creg] <- as.character(ridx)
}
outmat <- labsep(reglabels,reglist)
regmat <- outmat[[1]]
regsum <- outmat[[2]]
regcolors <- c("#E6F5C9","#0000FF","#AFEEEE","#8B4513","#00FFFF","#80B1D3","#2F4F4F",
               "#E31A1C","#377EB8","#FCCDE5","#A6D854","#9400D3","#808000","#FFFF33",
               "#33A02C","#FDB462","#FF7F00","#CCEBC5","#696969","#F0027F","#F5F5F5",
               "#DAA520")
names(regcolors) <- reglist

#Set parameters.
nroi <- 360
nk <- as.numeric(k)
states <- c('sc','sFC',paste0('s',1:nk))
nstates <- length(states)
ctrl_list <- c('ave','mod','deg')
nctrl <- length(ctrl_list)

# Read --------------------------------------------------------------------

#Read in the full Gramian error table to find states with errors.
infile <- paste0(basepath,'/state_images/SC_sFC_dFC_gram.csv')
gramall <- read.csv(infile,header=F)
gramstates <- gramall[(gramall[,2]!=0),1]
ngram <- length(gramstates)

#Go through each control type for values and limits.
ctrl_collect <- list()
maxmat <- matrix(NA,nctrl,nstates)
rownames(maxmat) <- ctrl_list
colnames(maxmat) <- states
minmat <- maxmat
for (cctrl in ctrl_list) {
  
  #Get SC.
  if (cctrl == 'deg') {
    infile <- paste0(scpath,'deg_sc.csv')
  } else {
    infile <- paste0(scpath,'orig_',cctrl,'_sc.csv')
  }
  predlabs <- paste0('sc_r',1:nroi)
  scmat <- read.csv(infile,row.names=1,header=F)
  colnames(scmat) <- predlabs
  
  #Remove SC errors.
  errstate <- 'SC'
  if (errstate %in% gramstates) {
    
    #Read regions to remove.
    infile <- paste0(basepath,'/state_images/',errstate,'_err.csv')
    errmat <- read.csv(infile,header=T)
    regrem <- paste0('sc_r',errmat[,'Region'])
    keeplabs <- predlabs[!(predlabs %in% regrem)]
    
    #Turn columns to NA.
    scmat[,regrem] <- NA
  }

  #Get sFC.
  if (cctrl == 'deg') {
    controltype <- 'abs_deg'
  } else {
    controltype <- cctrl
  }
  predlabs <- paste0('sFC_r',1:nroi)
  infile <- paste0(sFCpath,controltype,'_sFC.csv')
  sFCmat <- read.csv(infile,row.names=1,header=F)
  colnames(sFCmat) <- predlabs
  
  #Remove sFC errors.
  errstate <- 'sFC'
  if (errstate %in% gramstates) {
    
    #Read regions to remove.
    infile <- paste0(basepath,'/state_images/',errstate,'_err.csv')
    errmat <- read.csv(infile,header=T)
    regrem <- paste0('sFC_r',errmat[,'Region'])
    keeplabs <- predlabs[!(predlabs %in% regrem)]
    
    #Turn columns to NA.
    sFCmat[,regrem] <- NA
  }
  
  #Get dFC.
  if (cctrl == 'deg') {
    controltype <- 'abs_deg'
  } else {
    controltype <- cctrl
  }
  predlabs <- c()
  for (kidx in 1:nk) {
    predlabs <- c(predlabs,paste0('s',kidx,'_r',1:nroi))
  }
  infile <- paste0(basepath,controltype,'_tab.csv')
  dFCmat <- read.csv(infile,row.names=1,header=F)
  colnames(dFCmat) <- predlabs
  
  #Remove dFC errors.
  for (kidx in 1:nk) {
    errstate <- paste0('s',as.character(kidx))
    if (errstate %in% gramstates) {
      
      #Read regions to remove.
      infile <- paste0(basepath,'/state_images/',errstate,'_err.csv')
      errmat <- read.csv(infile,header=T)
      regrem <- paste0(errstate,'_r',errmat[,'Region'])
      keeplabs <- predlabs[!(predlabs %in% regrem)]
      
      #Turn columns to NA.
      dFCmat[,regrem] <- NA
    }
  }
  
  #Merge.
  predmat <- cbind(scmat,sFCmat,dFCmat)
  plotcollect <- matrix(NA,nrow=nroi,ncol=nstates)
  colnames(plotcollect) <- states
  for (cstate in states) {
    cmat <- predmat[,grep(cstate,colnames(predmat))]
    outkey <- paste0(cctrl,'.',cstate)
    cmean <- colMeans(cmat)
    names(cmean) <- NULL
    ctrl_collect[[outkey]] <- cmean
    maxmat[cctrl,cstate] <- max(cmean,na.rm=T)
    minmat[cctrl,cstate] <- min(cmean,na.rm=T)
  }
}

#Collapse dFC for limits.
setlist <- c('SC','sFC','dFC')
maxlim <- matrix(NA,nrow=3,ncol=3)
rownames(maxlim) <- ctrl_list
colnames(maxlim) <- setlist
minlim <- maxlim
maxlim[1:3,1:2] <- maxmat[1:3,1:2]
maxlim[1:3,3] <- apply(maxmat[1:3,3:nstates],1,max,na.rm=T)
minlim[1:3,1:2] <- minmat[1:3,1:2]
minlim[1:3,3] <- apply(minmat[1:3,3:nstates],1,min,na.rm=T)

# Plot --------------------------------------------------------------------

#Set parameters.
xwidth <- 3*3
xheight <- 2.2*(2+nk)
lheight <- 2.2*3
titlesize <- 15
wantcolor <- c('red','blue','green')

#Collect plots.
plist <- list()
pidx <- 1
for (cstate in states) {
  
  #Produce label.
  if (cstate=='sc') {
    cset <- 'SC'
    cstate_lab <- 'SC'
  } else if (cstate=='sFC') {
    cset <- 'sFC'
    cstate_lab <- 'sFC'
  } else {
    cset <- 'dFC'
    cstate_lab <- cstate
  }
  
  #Plot control values.
  for (cctrl in ctrl_list) {
 
    #Produce color.
    if (cctrl=='ave') {
      ccolor <- 'Reds'
    } else if (cctrl=='mod') {
      ccolor <- 'Blues'
    } else if (cctrl=='deg') {
      ccolor <- 'Greens'
    }
    
    #Extract.
    inkey <- paste0(cctrl,'.',cstate)
    plotmat <- ctrl_collect[[inkey]]
    vmat <- as.matrix(plotmat)
    colnames(vmat) <- c('Score')
    maxscore <- maxlim[cctrl,cset]
    minscore <- minlim[cctrl,cset]
    
    #Plot.
    curr_v <- cbind(atlasord,vmat)
    atlas_v <- inner_join(base_atlas,curr_v)
    ggobj <- atlas_v %>%
      ggplot() +
      geom_brain(atlas=glasser,
                 position=position_brain(hemi~side),
                 aes(fill=Score)) +
      scale_fill_distiller(palette=ccolor,
                        limits=c(minscore,maxscore),
                        na.value='grey') +
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
            plot.margin=grid::unit(c(1,0,1,0),'mm'))
    plist[[pidx]] <- ggobj
    pidx <- pidx + 1
  }
}

#Save.
outfile <- paste0(outpath,'threectrl.jpg')
ggsave(outfile,arrangeGrob(grobs=plist,ncol=3),
       units='in',width=xwidth,height=xheight)

#Save color bars for each.
llist <- list()
lidx <- 1
for (cset in setlist) {
  for (cctrl in ctrl_list) {
    
    #Produce color.
    if (cctrl=='ave') {
      ccolor <- 'Reds'
    } else if (cctrl=='mod') {
      ccolor <- 'Blues'
    } else if (cctrl=='deg') {
      ccolor <- 'Greens'
    }
    
    #Plot.
    maxscore <- maxlim[cctrl,cset]
    minscore <- minlim[cctrl,cset]
    vmat <- matrix(0,360,1)
    colnames(vmat) <- c('Score')
    curr_v <- cbind(atlasord,vmat)
    atlas_v <- inner_join(base_atlas,curr_v)
    ggobj <- atlas_v %>%
      ggplot() +
      geom_brain(atlas=glasser,
                 position=position_brain(hemi~side),
                 aes(fill=Score)) +
      scale_fill_distiller(palette=ccolor,
                           limits=c(minscore,maxscore),
                           na.value='grey') +
      theme_classic() +
      theme(plot.title = element_text(hjust=0.05),
            axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            title=element_blank(),
            legend.margin = margin(t=0,
                                   r=17,
                                   b=0,
                                   l=0, unit = "pt"),
            plot.margin=grid::unit(c(1,0,1,0),'mm'))
    llist[[lidx]] <- ggobj
    lidx <- lidx + 1
  }
}

#Save.
outfile <- paste0(outpath,'threectrl_legend.jpg')
ggsave(outfile,arrangeGrob(grobs=llist,ncol=3),
       units='in',width=xwidth,height=lheight)


# Detailed - State-wise ---------------------------------------------------

#Set parameters.
xwidth <- 3*6
xheight <- 2.2*(2+nk)
wantcolor <- c('red','blue','green')
nacolor <- '#222222'
qlist <- c(80,85,90,95,99)
nq <- length(qlist)

#Plot control values.
for (cq in qlist) {
  for (cctrl in ctrl_list) {
    
    #Produce color.
    if (cctrl=='ave') {
      ccolor <- 'Reds'
      radcolor <- 'red'
    } else if (cctrl=='mod') {
      ccolor <- 'Blues'
      radcolor <- 'red'
    } else if (cctrl=='deg') {
      ccolor <- 'Greens'
      radcolor <- 'red'
    }
  
    #Collect plots.
    plist <- list()
    pltidx <- 1
    for (cstate in states) {
      
      #Produce label.
      if (cstate=='sc') {
        cset <- 'SC'
        cstate_lab <- 'SC'
      } else if (cstate=='sFC') {
        cset <- 'sFC'
        cstate_lab <- 'sFC'
      } else {
        cset <- 'dFC'
        cstate_lab <- cstate
      }
      
      #Extract.
      inkey <- paste0(cctrl,'.',cstate)
      cmat <- ctrl_collect[[inkey]]
      qtile <- quantile(cmat,probs=(cq/100),na.rm=T)
      plotmat <- cmat
      plotmat[cmat<qtile] <- NA

      #Plot raw.
      vmat <- as.matrix(plotmat)
      colnames(vmat) <- c('Score')
      maxscore <- maxlim[cctrl,cset]
      minscore <- minlim[cctrl,cset]
      curr_v <- cbind(atlasord,vmat)
      atlas_v <- inner_join(base_atlas,curr_v)
      ggobj <- atlas_v %>%
        ggplot() +
        geom_brain(atlas=glasser,
                   position=position_brain(hemi~side),
                   aes(fill=Score)) +
        scale_fill_distiller(palette=ccolor,
                             limits=c(minscore,maxscore),
                             na.value=nacolor) +
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
              plot.margin=grid::unit(c(1,0,1,0),'mm'))
      plist[[pltidx]] <- ggobj
      pltidx <- pltidx + 1
      
      #Plot absolute.
      vmat <- as.matrix(plotmat)
      vmat[vmat>0] <- 0.25
      colnames(vmat) <- c('Score')
      curr_v <- cbind(atlasord,vmat)
      atlas_v <- inner_join(base_atlas,curr_v)
      ggobj <- atlas_v %>%
        ggplot() +
        geom_brain(atlas=glasser,
                   position=position_brain(hemi~side),
                   aes(fill=Score)) +
        scale_fill_distiller(palette=ccolor,
                             limits=c(0,1),
                             na.value=nacolor) +
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
              plot.margin=grid::unit(c(1,0,1,0),'mm'))
      plist[[pltidx]] <- ggobj
      pltidx <- pltidx + 1
      
      #Plot network plot.
      arealabels <- netlabels
      arealist <- netlist
      nareas <- nnet
      areamat <- netmat
      areasum <- netsum
      areacolors <- netcolors
      rad_txtsize <- 3
      outlist <- uni_areaplot(plotmat,cstate_lab,arealabels,arealist,nareas,areamat,areasum,areacolors,
                              rad_txtsize,atlasord,base_atlas,nacolor,radcolor,plist,pltidx) 
      plist <- outlist[[1]]
      pltidx <- outlist[[2]]
      
      #Plot region plot.
      arealabels <- reglabels
      arealist <- reglist
      nareas <- nreg
      areamat <- regmat
      areasum <- regsum
      areacolors <- regcolors
      rad_txtsize <- 2.5
      outlist <- uni_areaplot(plotmat,cstate_lab,arealabels,arealist,nareas,areamat,areasum,areacolors,
                              rad_txtsize,atlasord,base_atlas,nacolor,radcolor,plist,pltidx) 
      plist <- outlist[[1]]
      pltidx <- outlist[[2]]
    }
    
    #Save.
    outfile <- paste0(outpath,cctrl,'_',as.character(cq),'.jpg')
    ggsave(outfile,arrangeGrob(grobs=plist,ncol=6),
           units='in',width=xwidth,height=xheight)
  }
}
