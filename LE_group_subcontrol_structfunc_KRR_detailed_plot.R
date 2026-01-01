# For a given k for clustering, type of SC normalization, type and percentage 
# of thresholding, controllability type, the set of matrices, the set of
# cognitive variables, CV repetitions, CV inner fold number, CV outer fold number,
# and permutation number, plot the regions with significant Haufe scores and 
# corresponding radar plots for the networks and regions to analyze.
# Output:
# base_outfile,'.jpg' Brain plots with radar plots for networks and regions, where the radar is maximized.
# base_outfile,'_unmax.jpg' Brain plots with radar plots for networks and regions, where the radar is unchanged.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggseg)
library(ggsegGlasser)
library(scico)
library(RColorBrewer)
library(ggradar)

# RB Functions ------------------------------------------------------------

#Create red-blue brain plotter.
rbplot <- function(cstate,minscore,maxscore,atlasord,base_atlas,plist,pltidx) {
  
  #Plot.
  Score <- cstate
  curr_v <- cbind(atlasord,Score)
  atlas_v <- inner_join(base_atlas,curr_v)
  ggobj <- atlas_v %>%
    ggplot() +
    geom_brain(atlas=glasser, 
               position=position_brain(hemi~side),
               aes(fill=Score)) +
    scale_fill_scico(palette='berlin',
                     midpoint=0,
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
  plist[[pltidx]] <- ggobj
  pltidx <- pltidx + 1
  
  #Return items.
  return(list(plist,pltidx))
}

#Create red-blue brain plotter with extreme values.
rbXplot <- function(cstate,atlasord,base_atlas,rbcolors,plist,pltidx) {
  
  #Plot.
  Score <- cstate
  Score[Score > 0] <- 1
  Score[Score < 0] <- -1
  curr_v <- cbind(atlasord,Score)
  atlas_v <- inner_join(base_atlas,curr_v)
  ggobj <- atlas_v %>%
    ggplot() +
    geom_brain(atlas=glasser,
               position=position_brain(hemi~side),
               aes(fill=Score)) +
    scale_fill_gradientn(colors=rbcolors,
                         limits=c(-1,1),
                         breaks=c(-1,0,1),
                         labels=format(c(-1,0,1)),
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
  plist[[pltidx]] <- ggobj
  pltidx <- pltidx + 1
  
  #Return items.
  return(list(plist,pltidx))
}

#Generate red-blue legend.
rblegend <- function(minscore,maxscore,atlasord,base_atlas) {
  
  #Plot.
  Score <- c(rep(minscore,180),rep(maxscore,180))
  curr_v <- cbind(atlasord,Score)
  atlas_v <- inner_join(base_atlas,curr_v)
  ggobj <- atlas_v %>%
    ggplot() +
    geom_brain(atlas=glasser, 
               position=position_brain(hemi~side),
               aes(fill=Score)) +
    scale_fill_scico(palette='berlin',
                     midpoint=0,
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
          legend.margin = margin(t=0,
                                 r=17, 
                                 b=0, 
                                 l=0, unit = "pt"),
          plot.margin=grid::unit(c(1,0,1,0),'mm')) 
  
  #Return items.
  return(ggobj)
}

# Area Functions ----------------------------------------------------------

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

#Create area plotter with radar maximization.
areaplot <- function(cstate,cstatelab,arealabels,arealist,nareas,areamat,areasum,areacolors,
                     rad_txtsize,atlasord,base_atlas,nacolor,bcolor,plist,pltidx) {
  
  #Extract positive loadings for areas.
  cstate_plot <- cstate
  cstate_plot[cstate>0] <- 1
  cstate_plot[cstate<=0] <- NA
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
  
  #Extract negative loadings.
  cstate_plot <- cstate
  cstate_plot[cstate>=0] <- NA
  cstate_plot[cstate<0] <- 1
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
  
  #Extract positive and negative.
  pos.state <- cstate
  pos.state[cstate < 0] <- NA
  neg.state <- cstate
  neg.state[cstate > 0] <- NA
  neg.state <- -neg.state
  
  #Project to networks.
  pos.state.rep <- matrix(pos.state,nrow=360,ncol=nareas)
  pos.state.values <- pos.state.rep
  pos.state.values[!areamat] <- NA
  neg.state.rep <- matrix(neg.state,nrow=360,ncol=nareas)
  neg.state.values <- neg.state.rep
  neg.state.values[!areamat] <- NA
  
  #Find dice using count of regions of interest, count of regions in the network,
  #and count of the overlap.
  posA <- sum(!is.na(pos.state))
  posB <- areasum
  posAB <-  colSums(!is.na(pos.state.values))
  pos.state.dice <- dice(posAB,posA,posB)
  negA <- sum(!is.na(neg.state))
  negB <- areasum
  negAB <-  colSums(!is.na(neg.state.values))
  neg.state.dice <- dice(negAB,negA,negB)
  
  #Maximize nonzero positive and negative and plot.
  statepos <- sum((cstate>0),na.rm=T)
  stateneg <- sum((cstate<0),na.rm=T)
  ctitle <- paste0(toupper(cstatelab),' | P=',as.character(statepos),' N=',as.character(stateneg))
  pos.state.sep <- pos.state.dice
  neg.state.sep <-  neg.state.dice
  cposmax <- max(pos.state.sep)
  cnegmax <- max(neg.state.sep)
  if (cposmax > cnegmax) {
    cdiff <- cposmax - cnegmax
    neg.state.sep[neg.state.sep>0] <- neg.state.sep[neg.state.sep>0] + cdiff
  } else if (cposmax < cnegmax) {
    cdiff <- cnegmax - cposmax
    pos.state.sep[pos.state.sep>0] <- pos.state.sep[pos.state.sep>0] + cdiff
  }
  badder = data.frame('Group'=c('P','N'))
  both.sep <- rbind(pos.state.sep,neg.state.sep)
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
                   group.colours=bcolor,
                   legend.position='none') +
    theme(plot.margin=grid::unit(c(1,0,1,0),'mm'),
          plot.title=element_text(size=titlesize,hjust=0,face='bold')) +
    labs(title=ctitle)
  plist[[pltidx]] <- ggobj
  pltidx <- pltidx + 1
  
  #Find averages. Don't divide by zero, values with 1 will be divided by 1 to 
  # produce the value, values with 0 will be divided by 1 to produce 0.
  npos <- colSums(!is.na(pos.state.values))
  npos[npos==0] <- 1
  pos.state.avg <- t(data.frame(colSums(pos.state.values,na.rm=T)/npos))
  nneg <- colSums(!is.na(neg.state.values))
  nneg[nneg==0] <- 1
  neg.state.avg <- t(data.frame(colSums(neg.state.values,na.rm=T)/nneg))
  
  #Maximize nonzero positive and negative and plot.
  statepos <- round(mean(pos.state,na.rm=T),3)
  stateneg <- round(mean(neg.state,na.rm=T),3)
  if (is.na(statepos)) {statepos <- 0}
  if (is.na(stateneg)) {stateneg <- 0}
  ctitle <- paste0(toupper(cstatelab),' | P=',as.character(statepos),' N=',as.character(stateneg))
  pos.state.sep <- pos.state.avg
  neg.state.sep <-  neg.state.avg
  cposmax <- max(pos.state.sep)
  cnegmax <- max(neg.state.sep)
  if (cposmax > cnegmax) {
    cdiff <- cposmax - cnegmax
    neg.state.sep[neg.state.sep>0] <- neg.state.sep[neg.state.sep>0] + cdiff
  } else if (cposmax < cnegmax) {
    cdiff <- cnegmax - cposmax
    pos.state.sep[pos.state.sep>0] <- pos.state.sep[pos.state.sep>0] + cdiff
  }
  badder = data.frame('Group'=c('P','N'))
  both.sep <- rbind(pos.state.sep,neg.state.sep)
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
                   group.colours=bcolor,
                   legend.position='none') +
    theme(plot.margin=grid::unit(c(1,0,1,0),'mm'),
          plot.title=element_text(size=titlesize,hjust=0,face='bold')) +
    labs(title=ctitle)
  plist[[pltidx]] <- ggobj
  pltidx <- pltidx + 1
  
  #Return items.
  return(list(plist,pltidx))
}

#Get max-min for radar plots.
maxmin_areaplot <- function(cstate,arealabels,nareas,areamat,areasum) {
  
  #Extract positive and negative.
  pos.state <- cstate
  pos.state[cstate < 0] <- NA
  neg.state <- cstate
  neg.state[cstate > 0] <- NA
  neg.state <- -neg.state
  
  #Project to networks.
  pos.state.rep <- matrix(pos.state,nrow=360,ncol=nareas)
  pos.state.values <- pos.state.rep
  pos.state.values[!areamat] <- NA
  neg.state.rep <- matrix(neg.state,nrow=360,ncol=nareas)
  neg.state.values <- neg.state.rep
  neg.state.values[!areamat] <- NA
  
  #Find dice using count of regions of interest, count of regions in the network,
  #and count of the overlap.
  posA <- sum(!is.na(pos.state))
  posB <- areasum
  posAB <-  colSums(!is.na(pos.state.values))
  pos.state.dice <- dice(posAB,posA,posB)
  negA <- sum(!is.na(neg.state))
  negB <- areasum
  negAB <-  colSums(!is.na(neg.state.values))
  neg.state.dice <- dice(negAB,negA,negB)
  
  #Find averages. Don't divide by zero, values with 1 will be divided by 1 to 
  # produce the value, values with 0 will be divided by 1 to produce 0.
  npos <- colSums(!is.na(pos.state.values))
  npos[npos==0] <- 1
  pos.state.avg <- t(data.frame(colSums(pos.state.values,na.rm=T)/npos))
  nneg <- colSums(!is.na(neg.state.values))
  nneg[nneg==0] <- 1
  neg.state.avg <- t(data.frame(colSums(neg.state.values,na.rm=T)/nneg))
  
  #Combine.
  outmat <- rbind(pos.state.dice,neg.state.dice,pos.state.avg,neg.state.avg)
  rownames(outmat) <- c('pos_dice','neg_dice','pos_avg','neg_avg')
  
  #Return items.
  return(outmat)
}

#Create area plotter with no radar maximization.
unmax_areaplot <- function(cstate,cstatelab,arealabels,arealist,nareas,areamat,areasum,areacolors,
                           rad_txtsize,atlasord,base_atlas,nacolor,bcolor,plist,pltidx,
                           areadice,areaavg) {
  
  #Extract positive loadings for areas.
  cstate_plot <- cstate
  cstate_plot[cstate>0] <- 1
  cstate_plot[cstate<=0] <- NA
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
  
  #Extract negative loadings.
  cstate_plot <- cstate
  cstate_plot[cstate>=0] <- NA
  cstate_plot[cstate<0] <- 1
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
  
  #Extract positive and negative.
  pos.state <- cstate
  pos.state[cstate < 0] <- NA
  neg.state <- cstate
  neg.state[cstate > 0] <- NA
  neg.state <- -neg.state
  
  #Project to networks.
  pos.state.rep <- matrix(pos.state,nrow=360,ncol=nareas)
  pos.state.values <- pos.state.rep
  pos.state.values[!areamat] <- NA
  neg.state.rep <- matrix(neg.state,nrow=360,ncol=nareas)
  neg.state.values <- neg.state.rep
  neg.state.values[!areamat] <- NA
  
  #Find dice.
  posA <- sum(!is.na(pos.state))
  posB <- areasum
  posAB <-  colSums(!is.na(pos.state.values))
  pos.state.dice <- dice(posAB,posA,posB)
  negA <- sum(!is.na(neg.state))
  negB <- areasum
  negAB <-  colSums(!is.na(neg.state.values))
  neg.state.dice <- dice(negAB,negA,negB)
  
  #Put in obtained max and plot.
  statepos <- sum((cstate>0),na.rm=T)
  stateneg <- sum((cstate<0),na.rm=T)
  ctitle <- paste0(toupper(cstatelab),' | P=',as.character(statepos),' N=',as.character(stateneg))
  pos.state.sep <- pos.state.dice
  neg.state.sep <-  neg.state.dice
  badder = data.frame('Group'=c('P','N'))
  both.sep <- rbind(pos.state.sep,neg.state.sep)
  colnames(both.sep) <- colnames(areamat)
  both.sep <- cbind(badder,both.sep)
  maxrad <- areadice
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
                   group.colours=bcolor,
                   legend.position='none') +
    theme(plot.margin=grid::unit(c(1,0,1,0),'mm'),
          plot.title=element_text(size=titlesize,hjust=0,face='bold')) +
    labs(title=ctitle)
  plist[[pltidx]] <- ggobj
  pltidx <- pltidx + 1
  
  #Find averages.
  npos <- colSums(!is.na(pos.state.values))
  npos[npos==0] <- 1
  pos.state.avg <- t(data.frame(colSums(pos.state.values,na.rm=T)/npos))
  nneg <- colSums(!is.na(neg.state.values))
  nneg[nneg==0] <- 1
  neg.state.avg <- t(data.frame(colSums(neg.state.values,na.rm=T)/nneg))
  
  #Put in obtained max and plot.
  statepos <- round(mean(pos.state,na.rm=T),3)
  stateneg <- round(mean(neg.state,na.rm=T),3)
  if (is.na(statepos)) {statepos <- 0}
  if (is.na(stateneg)) {stateneg <- 0}
  ctitle <- paste0(toupper(cstatelab),' | P=',as.character(statepos),' N=',as.character(stateneg))
  pos.state.sep <- pos.state.avg
  neg.state.sep <-  neg.state.avg
  badder = data.frame('Group'=c('P','N'))
  both.sep <- rbind(pos.state.sep,neg.state.sep)
  colnames(both.sep) <- colnames(areamat)
  both.sep <- cbind(badder,both.sep)
  maxrad <- areaavg
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
                   group.colours=bcolor,
                   legend.position='none') +
    theme(plot.margin=grid::unit(c(1,0,1,0),'mm'),
          plot.title=element_text(size=titlesize,hjust=0,face='bold')) +
    labs(title=ctitle)
  plist[[pltidx]] <- ggobj
  pltidx <- pltidx + 1
  
  #Return items.
  return(list(plist,pltidx))
}

#Generate area plotter legend.
arealegend <- function(arealabels,arealist,areacolors,
                       atlasord,base_atlas) {
  
  #Extract loadings for areas.
  Score <- factor(arealist[as.numeric(arealabels)],levels=arealist)
  
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
          title=element_blank(),
          legend.margin = margin(t=0,
                                 r=17, 
                                 b=0, 
                                 l=0, unit = "pt"),
          plot.margin=grid::unit(c(1,0,1,0),'mm')) 
  
  #Return items.
  return(ggobj)
}

# Plotting Functions ------------------------------------------------------

#Put together the plotting functions for maximized radar plots.
max_plotter <- function(score_thres,minscore,maxscore,base_outfile,
                        atlasord,base_atlas,nacolor,bcolor,xwidth,xheight,
                        states,nstates,nroi,
                        rbcolors,
                        netlabels,netlist,nnet,netmat,netsum,netcolors,
                        reglabels,reglist,nreg,regmat,regsum,regcolors) {
  
  #Go through each state for maximized radar versions.
  plist <- list()
  pltidx <- 1
  for (stidx in 1:nstates) {
    
    #Extract.
    cstatelab <- states[stidx]
    cstatefeats <- paste0(cstatelab,'_r',1:nroi)
    cstate <- score_thres[cstatefeats]
    names(cstate) <- NULL
    
    #Get red-blue plot.
    outmat <- rbplot(cstate,minscore,maxscore,atlasord,base_atlas,plist,pltidx)
    plist <- outmat[[1]]
    pltidx <- outmat[[2]]
    
    #Get red-blue extreme plot.
    outmat <- rbXplot(cstate,atlasord,base_atlas,rbcolors,plist,pltidx)
    plist <- outmat[[1]]
    pltidx <- outmat[[2]]
    
    #Get area plots for networks.
    arealabels <- netlabels
    arealist <- netlist
    nareas <- nnet
    areamat <- netmat
    areasum <- netsum
    areacolors <- netcolors
    rad_txtsize <- 3
    outmat <- areaplot(cstate,cstatelab,arealabels,arealist,nareas,areamat,areasum,areacolors,
                       rad_txtsize,atlasord,base_atlas,nacolor,bcolor,plist,pltidx)
    plist <- outmat[[1]]
    pltidx <- outmat[[2]]
    
    #Get area plots for regions.
    arealabels <- reglabels
    arealist <- reglist
    nareas <- nreg
    areamat <- regmat
    areasum <- regsum
    areacolors <- regcolors
    rad_txtsize <- 2.5
    outmat <- areaplot(cstate,cstatelab,arealabels,arealist,nareas,areamat,areasum,areacolors,
                       rad_txtsize,atlasord,base_atlas,nacolor,bcolor,plist,pltidx)
    plist <- outmat[[1]]
    pltidx <- outmat[[2]]
  }
  
  #Save.
  outfile <- paste0(base_outfile,'.jpg')
  ggsave(outfile,arrangeGrob(grobs=plist,ncol=10),
         units='in',width=xwidth,height=xheight)
}

#Put together the plotting functions for unmaximized radar plots.
unmax_plotter <- function(score_thres,minscore,maxscore,base_outfile,
                          atlasord,base_atlas,nacolor,bcolor,xwidth,xheight,
                          states,nstates,nroi,
                          rbcolors,
                          netlabels,netlist,nnet,netmat,netsum,netcolors,
                          reglabels,reglist,nreg,regmat,regsum,regcolors) {
  
  #Go through each state to retrieve max-min values.
  netmax <- matrix(NA,4,nstates)
  rownames(netmax) <- c('pos_dice','neg_dice','pos_avg','neg_avg')
  colnames(netmax) <- states
  regmax <- netmax
  for (stidx in 1:nstates) {
    
    #Extract.
    cstatelab <- states[stidx]
    cstatefeats <- paste0(cstatelab,'_r',1:nroi)
    cstate <- score_thres[cstatefeats]
    names(cstate) <- NULL
    
    #Get table of network radar values and insert.
    arealabels <- netlabels
    nareas <- nnet
    areamat <- netmat
    areasum <- netsum
    outrad <- maxmin_areaplot(cstate,arealabels,nareas,areamat,areasum)
    netmax[,cstatelab] <- do.call('pmax',data.frame(outrad))
    
    #Get table of region radar values and insert.
    arealabels <- reglabels
    nareas <- nreg
    areamat <- regmat
    areasum <- regsum
    outrad <- maxmin_areaplot(cstate,arealabels,nareas,areamat,areasum)
    regmax[,cstatelab] <- do.call('pmax',data.frame(outrad))
  }
  
  #Get the largest value for each metric.
  netdice <- max(netmax[1:2,])
  netavg <- max(netmax[3:4,])
  regdice <- max(regmax[1:2,])
  regavg <- max(regmax[3:4,])
  
  #Go through each state for unmaximized radar versions.
  plist <- list()
  pltidx <- 1
  for (stidx in 1:nstates) {
    
    #Extract.
    cstatelab <- states[stidx]
    cstatefeats <- paste0(cstatelab,'_r',1:nroi)
    cstate <- score_thres[cstatefeats]
    names(cstate) <- NULL
    
    #Get red-blue plot.
    outmat <- rbplot(cstate,minscore,maxscore,atlasord,base_atlas,plist,pltidx)
    # outmat <- rbplot(cstate,minscore,maxscore,atlasord,base_atlas,nacolor,rbcolors,plist,pltidx)
    plist <- outmat[[1]]
    pltidx <- outmat[[2]]
    
    #Get red-blue extreme plot.
    outmat <- rbXplot(cstate,atlasord,base_atlas,rbcolors,plist,pltidx)
    # outmat <- rbXplot(cstate,atlasord,base_atlas,nacolor,rbcolors,plist,pltidx)
    plist <- outmat[[1]]
    pltidx <- outmat[[2]]
    
    #Get area plots for networks.
    arealabels <- netlabels
    arealist <- netlist
    nareas <- nnet
    areamat <- netmat
    areasum <- netsum
    areacolors <- netcolors
    rad_txtsize <- 3
    areadice <- netdice
    areaavg <- netavg
    outmat <- unmax_areaplot(cstate,cstatelab,arealabels,arealist,nareas,areamat,areasum,areacolors,
                             rad_txtsize,atlasord,base_atlas,nacolor,bcolor,plist,pltidx,
                             areadice,areaavg)
    plist <- outmat[[1]]
    pltidx <- outmat[[2]]
    
    #Get area plots for regions.
    arealabels <- reglabels
    arealist <- reglist
    nareas <- nreg
    areamat <- regmat
    areasum <- regsum
    areacolors <- regcolors
    rad_txtsize <- 2.5
    areadice <- regdice
    areaavg <- regavg
    outmat <- unmax_areaplot(cstate,cstatelab,arealabels,arealist,nareas,areamat,areasum,areacolors,
                             rad_txtsize,atlasord,base_atlas,nacolor,bcolor,plist,pltidx,
                             areadice,areaavg)
    plist <- outmat[[1]]
    pltidx <- outmat[[2]]
  }
  
  #Save.
  # do.call('grid.arrange',c(plist,ncol=10))
  outfile <- paste0(base_outfile,'_unmax.jpg')
  ggsave(outfile,arrangeGrob(grobs=plist,ncol=10),
         units='in',width=xwidth,height=xheight)
}


# Driver ------------------------------------------------------------------

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
nperm <- args[11]
print(paste(k,sctype,threstype,thresval,controltype,statetype,septype,nrep,
            inner_k,outer_k,nperm,sep=' '))

#Set base path and predict path.
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
                  subgroup,'/',k,'/SC_dFC/',sc_subgroup,'/collect/',threstype,'/',thresval,'/')
compath <- paste0(basepath,'KRRXFS/controlcomp_statecomp_cogcomp_',
                  nrep,'_',inner_k,'_',outer_k,'_',sctype,'/')
krrpath <- paste0(basepath,'KRRXFS/',controltype,'_',statetype,'_',septype,'_',
                  nrep,'_',inner_k,'_',outer_k,'_',sctype,'/')

#Set parameters.
nk <- as.numeric(k)
nroi <- 360

#Generate versions of interest for information.
full_featlabs <- c()
if (statetype == 'SC_sFC_dFCcat') {
  states <- c('sc','sFC',paste0('s',1:nk))
  full_featlabs <- c(full_featlabs,paste0('sc_r',1:nroi))
  full_featlabs <- c(full_featlabs,paste0('sFC_r',1:nroi))
  for (x in paste0('s',1:nk)) {full_featlabs <- c(full_featlabs,paste0(x,'_r',1:nroi))}
} else if (statetype == 'SC_dFCcat') {
  states <- c('sc',paste0('s',1:nk))
  full_featlabs <- c(full_featlabs,paste0('sc_r',1:nroi))
  for (x in paste0('s',1:nk)) {full_featlabs <- c(full_featlabs,paste0(x,'_r',1:nroi))}
} else if (statetype == 'dFCcat') {
  states <- paste0('s',1:nk)
  for (x in paste0('s',1:nk)) {full_featlabs <- c(full_featlabs,paste0(x,'_r',1:nroi))}
} else if (statetype == 'SC') {
  states <- c('sc')
  full_featlabs <- c(full_featlabs,paste0('sc_r',1:nroi))
} else if (statetype == 'sFC') {
  states <- c('sFC')
  full_featlabs <- c(full_featlabs,paste0('sFC_r',1:nroi))
}
nstates <- length(states)
nfull <- length(full_featlabs)
if (septype=='comCFAng') {
  coglist <- c('gCFA','P24_CR','PV','gFngCFA','gCngCFA')
  ncog <- length(coglist)
}

#Generate versions of interest for prediction.
qlist <- c(0,50,60,70,80,90)
nq <- length(qlist)

#Set plot parameters.
xwidth <- 3*10
xheight <- 0
inc <- 2.2
if (grepl('dFCcat',statetype,fixed=T)) {
  xheight <- xheight + (inc*nk)
}
if (grepl('SC',statetype,fixed=T)) {
  xheight <- xheight + inc
}
if (grepl('sFC',statetype,fixed=T)) {
  xheight <- xheight + inc
}
nacolor <- '#222222'
bcolor <- c('blue','red')
titlesize <- 15

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

#Set red and blue colors.
rbcolors <- c('#1f77b4','#222222','#d62728')

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

# Collect Data ------------------------------------------------------------
  

#Read in cross-cognition p-values.
infile <- paste0(compath,controltype,'_',statetype,'_',septype,'_TwoP_BHFDR_covha_feat.csv')
twop_BH <- read.csv(infile,row.names=1)

#Extract feature labels.
featlabs <- rownames(twop_BH)

#P versions.
for (pidx in 1:npvers) {
  cpver <- pvers[pidx]
  
  #Quantile versions.
  for (qidx in 1:nq) {
    cq <- qlist[qidx]
    
    #Produce column names for all versions of interest, including state and cognition.
    colcollect <- c()
    for (cstate in states) {
      for (ccog in coglist) {
        colcollect <- c(colcollect,paste0(cstate,'.',ccog))
      }
    }
    ncolcollect <- length(colcollect)
    
    #Collect values and maximum and minimum without thresholding.
    ccollector <- matrix(NA,nroi,ncolcollect)
    rownames(ccollector) <- paste0('r',1:nroi)
    colnames(ccollector) <- colcollect
    cmaxmin <- matrix(NA,2,ncolcollect)
    rownames(cmaxmin) <- c('max','min')
    colnames(cmaxmin) <- colcollect
    
    #Collect cognitive versions of interest.
    for (cidx in 1:ncog) {
      ccog <- coglist[cidx]
      
      #Read in the omnibus matrix.
      infile <- paste0(krrpath,ccog,'_covha_featscores.csv')
      omnimat <- read.csv(infile,row.names=1)
      
      #Select score to display.
      full_score <- omnimat[,'Score']
      
      #Select p-values.
      pval <- twop_BH[,cidx]
      
      #Threshold.
      score <- full_score
      score[pval >= 0.05] <- NA
      
      #Remove bottom X%.
      cq <- qlist[1]
      cqmat <- score[score>0]
      posq <- quantile(cqmat,probs=(cq/100),na.rm=T)
      cqmat <- -score[score<0]
      negq <- quantile(cqmat,probs=(cq/100),na.rm=T)
      score_thres <- score
      score_thres[(score>0)&(score<posq)] <- NA
      score_thres[(score<0)&(score>(-negq))] <- NA
      score <- score_thres
      
      #Add labels.
      names(full_score) <- rownames(omnimat)
      names(score) <- rownames(omnimat)
      
      #State types.
      for (cstate in states) {
        
        #Generate the collector label.
        collect_lab <- paste0(cstate,'.',ccog)
        
        #Generate region scores and labels.
        reg_score <- score[grepl(cstate,names(score))]
        wholenames <- names(reg_score)
        barenames <- gsub(paste0(cstate,'_'),'',wholenames)
        names(reg_score) <- barenames
        ccollector[barenames,collect_lab] <- reg_score
        
        #Generate max and min.
        reg_maxmin <- full_score[grepl(cstate,names(full_score))]
        cmaxmin['max',collect_lab] <- max(reg_maxmin,na.rm=T)
        cmaxmin['min',collect_lab] <- min(reg_maxmin,na.rm=T)
      }
    }
    
    #Collect values and max-min for comparisons.
    ccols <- c()
    for (ccog in coglist) {ccols <- c(ccols,ccog)}
    nccols <- length(ccols)
    threscollect <- matrix(NA,nfull,nccols)
    rownames(threscollect) <- full_featlabs
    colnames(threscollect) <- ccols
    thresmaxmin <- matrix(NA,2,nccols)
    rownames(thresmaxmin) <- c('max','min')
    colnames(thresmaxmin) <- ccols
    for (ccog in coglist) {
      ccol_lab <- ccog
      limvals <- c()
      for (cstate in states) {
        collect_lab <- paste0(cstate,'.',ccog)
        cvals <- ccollector[,collect_lab]
        clim <- cmaxmin[,collect_lab]
        wholenames <- paste0(cstate,'_r',1:nroi)
        threscollect[wholenames,ccol_lab] <- cvals
        limvals <- c(limvals,clim)
      }
      thresmaxmin['max',ccol_lab] <- max(limvals)
      thresmaxmin['min',ccol_lab] <- min(limvals)
    }
    
    #Generate path.
    plotpath <- paste0(krrpath,'plots_detailed/covha_BH_TwoP/')
    dir.create(plotpath,recursive=T)
    
    #For each cognitive version, plot corresponding versions.
    for (ccog in coglist) {
      
      #Raw feature.
      score_thres <- threscollect[,ccog]
      base_outfile <- paste0(plotpath,'covha_BH_TwoP_',ccog,'_',as.character(cq))
      minscore <- thresmaxmin['min',ccog]
      maxscore <- thresmaxmin['max',ccog]
      max_plotter(score_thres,minscore,maxscore,base_outfile,
                  atlasord,base_atlas,nacolor,bcolor,xwidth,xheight,
                  states,nstates,nroi,
                  rbcolors,
                  netlabels,netlist,nnet,netmat,netsum,netcolors,
                  reglabels,reglist,nreg,regmat,regsum,regcolors)
      unmax_plotter(score_thres,minscore,maxscore,base_outfile,
                    atlasord,base_atlas,nacolor,bcolor,xwidth,xheight,
                    states,nstates,nroi,
                    rbcolors,
                    netlabels,netlist,nnet,netmat,netsum,netcolors,
                    reglabels,reglist,nreg,regmat,regsum,regcolors)
      
      #Generate legend for the raw.
      ggobj <- rblegend(minscore,maxscore,atlasord,base_atlas) 
      outfile <- paste0(base_outfile,'_legend.jpg')
      ggsave(outfile,ggobj)
    }
  }
}

#Generate network and region legends.
plotpath <- paste0(krrpath,'plots_detailed/')
outfile <- paste0(plotpath,'ICN_legend.jpg')
ggobj <- arealegend(netlabels,netlist,netcolors,atlasord,base_atlas)
ggsave(outfile,ggobj,units='in',width=8,height=4.4)
outfile <- paste0(plotpath,'Region_legend.jpg')
ggobj <- arealegend(reglabels,reglist,regcolors,atlasord,base_atlas)
ggsave(outfile,ggobj,units='in',width=8,height=4.4)
