# For a given k for clustering, type of SC normalization, type and percentage 
# of thresholding, controllability type, the set of matrices, the set of
# cognitive variables, CV repetitions, CV inner fold number, CV outer fold number,
# and permutation number, plot the regions with significant FDR-corrected Haufe scores.
# Output:
# ccog,'_featscores_FDR.png' Brain plots for Haufe scores, thresholded for FDR p < 0.05.

# Set libraries.
library(tidyverse)
library(hdf5r)
library(ggplot2)
library(gridExtra)
# install.packages('sf')
library(ggseg)
library(ggh4x)
# install.packages('pak')
# pak::pak("ggsegverse/ggsegGlasser")
library(ggsegGlasser)
library(scico)

# Helper Functions --------------------------------------------------------

# Create label separator, 1 is within label.
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

# Driver ------------------------------------------------------------------

# Catches arguments.
args <- commandArgs(trailingOnly=T)
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
cat(k,sctype,threstype,thresval,controltype,statetype,septype,
    nrep,inner_k,outer_k,nperm,'\n')

# Set paths.
subgroup <- 'full'
sc_subgroup <- 'dr_full'
basepath <- paste0('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',
        subgroup , '/' , k ,
        '/SC_dFC/' , sc_subgroup , '/collect/' , threstype , '/' ,
        thresval)

# Set parameters.
nk <- as.numeric(k)
nroi <- 360
nconn <- nroi*(nroi-1)/2
ccolors <- c('#00A2FF','#FFFFFF','#FF0000')
  
# Read in our atlas organization.
atlasfile <- '/projects/jng/FLEXCOG/inputs/data/atlas/atlasorder.csv'
atlasord <- read_csv(atlasfile)[,'labelname'] %>%
  mutate(label=case_when(
    startsWith(labelname,'R_') ~ paste0('rh_',labelname),
    startsWith(labelname,'L_') ~ paste0('lh_',labelname)
  ))

# Plot --------------------------------------------------------------------

# Plot AC, MC, and S side by side with separate color bars for each for g, gF, and gC.
septype <- 'comCFAng'
coglist <- c('gCFA','P24_CR','PV')
ncog <- length(coglist)
ctrl_list <- c('ave','mod','abs_deg')
nctrl <- length(ctrl_list)
ctrl_labs <- c('Average Controllability','Modal Controllability','Strength')
statetype <- 'SC_dFCcat'
if (statetype=='SC_dFCcat') {
  statelist <- c('sc',paste0('s',as.character(1:nk)))
  nstates <- length(statelist)
}
for (cogidx in 1:ncog) {
  ccog <- coglist[cogidx]
  
  # Set dimensions.
  xinc <- 3
  xwidth <- xinc*3
  yheight <- 0
  yinc <- 2.13
  if (grepl('dFCcat',statetype,fixed=T)) {
    yheight <- yheight + (yinc*nk)
  }
  if (grepl('SC',statetype,fixed=T)) {
    yheight <- yheight + yinc
  }
  if (grepl('sFC',statetype,fixed=T)) {
    yheight <- yheight + yinc
  }
  yheight <- yheight + yinc
  
  # Collect plots.
  plist <- list()
  pidx <- 1
  lim_mat <- matrix(NA,nrow=nctrl,ncol=2)
  rownames(lim_mat) <- ctrl_list
  colnames(lim_mat) <- c('Min','Max')
  leglist <- list()
  for (ctrlidx in 1:nctrl) {
    cctrl <- ctrl_list[ctrlidx]
    cctrl_lab <- ctrl_labs[ctrlidx]
    
    # Read in the feature importance and p-values, then threshold.
    inpath <- paste0(basepath,'/KRRXFS/',cctrl,'_',statetype,'_',septype,'_',
                nrep , '_' , inner_k , '_' , outer_k , '_' , sctype)
    infile <- paste0(inpath,'/',ccog,'_covha_featscores_FDR.csv')
    inmat <- read.csv(infile,row.names=1)
    inmat <- inmat[inmat[,'TwoP_FDR']<0.05,]
    featlabs <- row.names(inmat)
    featall <- inmat[,'Score']
  
    # Get limits and append.
    clim <- c(min(featall),max(featall))
    lim_mat[cctrl,] <- clim
    
    # Do states.
    for (stidx in 1:nstates) {
      
      # Extract state.
      cstate <- statelist[stidx]
      stateyes <- grepl(cstate,featlabs,fixed=T)
      statefeat <- featall[stateyes]
      statelabs <- gsub(paste0(cstate,'_r'),'',featlabs[stateyes])
      names(statefeat) <- statelabs
      
      # Populate.
      vmat <- matrix(NA,nrow=nroi,ncol=1)
      rownames(vmat) <- as.character(c(1:nroi))
      vmat[statelabs,] <- statefeat
      colnames(vmat) <- 'Score'
      vplot <- atlasord %>%
        mutate(value=vmat[,'Score'])
      
      # Add plot and increment.
      border_sides <- if (stidx != nstates) 'b' else 'b'
      canno <- if (cctrl=='ave') toupper(cstate) else ''
      ctitle <- if (cstate=='sc') cctrl_lab else ''
      p <- ggplot() +
        geom_brain(atlas = glasser(),
                   mapping = aes(fill = value),
                   data = vplot,
                   view = c('lateral','medial'),
                   position = position_brain(view ~ hemi),
                   show.legend = F) +
        scale_fill_gradient2(low=ccolors[1],
                             mid=ccolors[2],
                             high=ccolors[3],
                             midpoint=0,
                             limits=clim,
                             na.value='black') +
        theme_brain() +
        theme(plot.title=element_text(hjust=0.5,vjust=-1.5,size=12,
                                      face='bold',family='Times New Roman'),
              plot.tag.position=c(0.09,0.91),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              plot.margin=grid::unit(c(0,0,0,0),'mm'),
              panel.border=element_blank(),
              panel.background = element_part_rect(side = border_sides,
                                                   color = 'black',
                                                   linewidth = 0.5,
                                                   fill = NA)) +
        ggtitle(ctitle) +
        annotate('text', x = -Inf, y = -Inf,
                 label = canno,
                 hjust = -0.1, vjust = -0.3,
                 size = 4.5, fontface = 'bold',family='Times New Roman')
      plist[[pidx]] <- p
      pidx <- pidx + 1
    }
    
    # Add legend.
    p <- ggplot() +
      geom_brain(atlas = glasser(),
                 mapping = aes(fill = value),
                 data = vplot,
                 view = c('lateral','medial'),
                 position = position_brain(view ~ hemi),
                 show.legend = T,
                 fill=NA,
                 colour=NA) +
      scale_fill_gradient2(low=ccolors[1],
                           mid=ccolors[2],
                           high=ccolors[3],
                           midpoint=0,
                           limits=clim) +
      theme_brain() +
      theme(legend.position      = 'top',
            legend.justification = 'center',
            legend.box           = 'vertical',
            legend.box.just      = 'right',
            legend.text          = element_text(size = 8,angle=90),
            legend.title         = element_text(size = 9,hjust=1,vjust=0.8,
                                                face='bold',family='Times New Roman'),
            legend.margin        = margin(t = 5, b = 5, unit = 'pt'),
            legend.spacing.y     = unit(-4, 'pt'),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            plot.margin       = grid::unit(c(0, 0, 0, 0), 'mm'))
    plist[[pidx]] <- p
    pidx <- pidx + 1
  }
    
  # Reorder, plot them all together, and save.
  plist_reordered <- vector('list', nctrl * nstates)
  for (ctrlidx in 1:nctrl) {
    for (stidx in 1:(nstates+1)) {
      new_idx <- (stidx - 1) * nctrl + ctrlidx
      old_idx <- (ctrlidx - 1) * (nstates+1) + stidx
      plist_reordered[[new_idx]] <- plist[[old_idx]]
    }
  }
  outpath <- paste0(basepath,'/KRRXFS/controlcomp_statecomp_cogcomp_',
                    nrep , '_' , inner_k , '_' , outer_k , '_' , sctype,'/plots')
  outfile <- paste0(outpath,'/',ccog,'_featscores_FDR.png')
  ggsave(outfile,arrangeGrob(grobs=plist_reordered,nrow=(nstates+1),ncol=nctrl,padding=unit(0,'mm')),
         units='in',width=xwidth,height=yheight,dpi=1440,bg='white')
}
print('Saved.')
