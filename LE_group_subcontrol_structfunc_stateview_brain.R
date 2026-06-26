# For the given k for clustering, type of SC normalization, and type and percentage of 
# thresholding, plot average controllability, modal controllability, and degree
# on the brain.
# Output:
# allbrain_controlplot.jpg Colored brain plots for average controllability, modal controllability, and degree for SC, sFC, and dFC.

#Set libraries.
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggseg)
library(ggh4x)
library(ggsegGlasser)
library(ggnewscale)
library(scico)

# Helper Functions --------------------------------------------------------

# Map state key to limit group.
get_lim_group <- function(cstate) {
  if (cstate == 'sc')  return('sc')
  if (cstate == 'sFC') return('sFC')
  return('dFC')
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

# Read in our atlas organization.
atlasfile <- '/projects/jng/FLEXCOG/inputs/data/atlas/atlasorder.csv'
atlasord <- read_csv(atlasfile)[,'labelname'] %>%
  mutate(label=case_when(
    startsWith(labelname,'R_') ~ paste0('rh_',labelname),
    startsWith(labelname,'L_') ~ paste0('lh_',labelname)
  ))

#Set parameters.
nroi <- 360
nk <- as.numeric(k)
statelist <- c('SC','sFC',paste0('s',as.character(1:nk)))
nstates <- length(statelist)

#Read in the full Gramian error table to find states with errors.
infile <- paste0(basepath,'SC_sFC_dFC_gram.csv')
gramall <- read.csv(infile,row.names=1)
gramstates <- rownames(gramall)
gramstates <- gramstates[gramall != 0]
ngram <- length(gramstates)

#Read SC AC, MC, and AD.
infile <- paste0(scpath,'orig_ave_sc.csv')
scave <- read.csv(infile,row.names=1,header=F)
infile <- paste0(scpath,'orig_mod_sc.csv')
scmod <- read.csv(infile,row.names=1,header=F)
infile <- paste0(scpath,'deg_sc.csv')
scdeg <- read.csv(infile,row.names=1,header=F)
predlabs <- paste0('sc_r',as.character(1:nroi))
colnames(scave) <- predlabs
colnames(scmod) <- predlabs
colnames(scdeg) <- predlabs

#If the error list contains SC.
errstate <- 'SC'
if (errstate %in% gramstates) {
  
  #Read in the regions to remove.
  infile <- paste0(basepath,'/state_images/',errstate,'_err.csv')
  errmat <- read.csv(infile)
  regrem <- errmat[,'Region']
  regrem <- paste0(errstate,'_r',as.character(regrem))
  scave[,regrem] <- NA
  scmod[,regrem] <- NA
}

#Average across subjects.
scave <- colMeans(scave)
scmod <- colMeans(scmod)
scdeg <- colMeans(scdeg)

#Read sFC AC, MC, and AD.
infile <- paste0(sFCpath,'ave_sFC.csv')
sFCave <- read.csv(infile,row.names=1,header=F)
infile <- paste0(sFCpath,'mod_sFC.csv')
sFCmod <- read.csv(infile,row.names=1,header=F)
infile <- paste0(sFCpath,'abs_deg_sFC.csv')
sFCdeg <- read.csv(infile,row.names=1,header=F)
predlabs <- paste0('sFC_r',as.character(1:nroi))
colnames(sFCave) <- predlabs
colnames(sFCmod) <- predlabs
colnames(sFCdeg) <- predlabs

#If the error list contains sFC.
errstate <- 'sFC'
if (errstate %in% gramstates) {
  
  #Read in the regions to remove.
  infile <- paste0(basepath,'/state_images/',errstate,'_err.csv')
  errmat <- read.csv(infile)
  regrem <- errmat[,'Region']
  regrem <- paste0(errstate,'_r',as.character(regrem))
  sFCave[,regrem] <- NA
  sFCmod[,regrem] <- NA
}

#Average across subjects.
sFCave <- colMeans(sFCave)
sFCmod <- colMeans(sFCmod)
sFCdeg <- colMeans(sFCdeg)

#Read dFC AC, MC, and AD.
infile <- paste0(basepath,'ave_tab.csv')
dFCave <- read.csv(infile,row.names=1,header=F)
infile <- paste0(basepath,'mod_tab.csv')
dFCmod <- read.csv(infile,row.names=1,header=F)
infile <- paste0(basepath,'abs_deg_tab.csv')
dFCdeg <- read.csv(infile,row.names=1,header=F)
predlabs <- c()
for (kidx in 1:nk) {
  for (ridx in 1:nroi) {
    cregion <- paste0('s',as.character(kidx),'_r',as.character(ridx))
    predlabs <- c(predlabs,cregion)
  }
}
colnames(dFCave) <- predlabs
colnames(dFCmod) <- predlabs
colnames(dFCdeg) <- predlabs

#If the error list contains dFC states.
for (kidx in 1:nk) {
  errstate <- paste0('s',as.character(kidx))
  if (errstate %in% gramstates) {
    
    #Read in the regions to remove.
    infile <- paste0(basepath,'/state_images/',errstate,'_err.csv')
    errmat <- read.csv(infile)
    regrem <- errmat[,'Region']
    regrem <- paste0(errstate,'_r',as.character(regrem))
    dFCave[,regrem] <- NA
    dFCmod[,regrem] <- NA
  }
}

#Average across subjects.
dFCave <- colMeans(dFCave)
dFCmod <- colMeans(dFCmod)
dFCdeg <- colMeans(dFCdeg)

# Plot --------------------------------------------------------------------

# Parameters.
ctrl_list <- c('ave','mod','deg')
ctrl_labs <- c('Average Controllability','Modal Controllability','Strength')
nctrl <- length(ctrl_list)
ctrl_colors <- c(ave = '#FF0000', mod = '#00A2FF', deg = '#097969')

# Build named list of data: each entry is a nroi-length named vector.
raw_data <- list()

# SC.
raw_data[['sc']] <- list(ave = scave, mod = scmod, deg = scdeg)

# sFC.
raw_data[['sFC']] <- list(ave = sFCave, mod = sFCmod, deg = sFCdeg)

# dFC states.
for (kidx in 1:nk) {
  cstate <- paste0('s', as.character(kidx))
  state_idx <- grepl(paste0('^', cstate, '_r'), names(dFCave))
  raw_data[[cstate]] <- list(
    ave = dFCave[state_idx],
    mod = dFCmod[state_idx],
    deg = dFCdeg[state_idx]
  )
}

# State list and layout: SC, sFC, then dFC states.
statelist_plot <- c('sc','sFC',paste0('s',as.character(1:nk)))
nstates_plot <- length(statelist_plot)

# Compute limits per control type. SC and sFC get their own limits; dFC states share limits.
lim_list <- list()
for (ctrlidx in 1:nctrl) {
  cctrl <- ctrl_list[ctrlidx]
  sc_vals  <- raw_data[['sc']][[cctrl]]
  sFC_vals <- raw_data[['sFC']][[cctrl]]
  dFC_vals <- unlist(lapply(paste0('s', 1:nk), function(s) raw_data[[s]][[cctrl]]))
  lim_list[[cctrl]] <- list(
    sc  = c(min(sc_vals,  na.rm = T), max(sc_vals,  na.rm = T)),
    sFC = c(min(sfc_vals, na.rm = T), max(sFC_vals, na.rm = T)),
    dFC = c(min(dfc_vals, na.rm = T), max(dFC_vals, na.rm = T))
  )
}

# Plot dimensions.
xinc <- 3
xwidth  <- xinc * nctrl
yinc    <- 2.13
yheight <- yinc * (2 + nk + 1)

# Build plot list, column-first then reorder.
plist <- list()
pidx  <- 1
for (ctrlidx in 1:nctrl) {
  cctrl <- ctrl_list[ctrlidx]
  cctrl_lab <- ctrl_labs[ctrlidx]
  ccolor <- ctrl_colors[cctrl]
  for (stidx in 1:nstates_plot) {
    cstate <- statelist_plot[stidx]
    lim_group <- get_lim_group(cstate)
    clim <- lim_list[[cctrl]][[lim_group]]
    
    # Extract values and strip state prefix from names.
    vals <- raw_data[[cstate]][[cctrl]]
    roi_idx <- gsub('.*_r', '', names(vals))
    
    # Populate.
    vmat <- matrix(NA, nrow = nroi, ncol = 1)
    rownames(vmat) <- as.character(1:nroi)
    vmat[roi_idx, ] <- vals
    colnames(vmat) <- 'Value'
    vplot <- atlasord %>%
      mutate(value = vmat[, 'Value'])
    border_sides <- if (stidx != nstates_plot) 'b' else 'b'
    cstate_lab <- if (cstate!='sFC') toupper(cstate) else 'sFC'
    canno <- if (cctrl=='ave') cstate_lab else ''
    if (stidx == 1) {ctitle <- cctrl_lab} else {ctitle <- ''}
    p <- ggplot() +
      geom_brain(atlas = glasser(),
                 mapping = aes(fill = value),
                 data = vplot,
                 view = c('lateral', 'medial'),
                 position = position_brain(view ~ hemi),
                 show.legend = F) +
      scale_fill_gradient(low = '#FFFFFF',
                          high = ccolor,
                          limits = clim,
                          na.value = 'black') +
      theme_brain() +
      ggtitle(ctitle) +
      theme(plot.title = element_text(hjust = 0.5, vjust = -1.5, size = 12,
                                      face='bold',family='Times New Roman'),
            plot.tag.position = c(0.09, 0.91),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            plot.margin = grid::unit(c(0, 0, 0, 0), 'mm'),
            panel.background = element_part_rect(side = border_sides,
                                               color = 'black',
                                               linewidth = 0.5,
                                               fill = NA,inherit.blank=F)) +
      annotate('text', x = -Inf, y = -Inf,
               label = canno,
               hjust = -0.1, vjust = -0.3,
               size = 4.5, fontface = 'bold',family='Times New Roman') 
    plist[[pidx]] <- p
    pidx <- pidx + 1
  }
  
  # Legend.
  dummy <- tibble(x = 1, y = 1:3, v = NA_real_)
  p_leg <- ggplot() +
    geom_tile(data    = dummy %>% filter(y == 1),
              mapping = aes(x = x, y = y, fill = v)) +
    scale_fill_gradient(low = '#FFFFFF',
                        high = ccolor,
                        limits = lim_list[[cctrl]][['sc']],
                        na.value = NA,
                        name = 'SC',
                        labels = function(x) sub("^0\\.", ".", as.character(x)),
                        guide = guide_colorbar(order = 1,
                                              direction = 'horizontal',
                                              title.position = 'left')) +
    new_scale_fill() +
    geom_tile(data = dummy %>% filter(y == 2),
              mapping = aes(x = x, y = y, fill = v)) +
    scale_fill_gradient(low = '#FFFFFF',
                        high = ccolor,
                        limits = lim_list[[cctrl]][['sFC']],
                        na.value = NA,
                        name = 'sFC',
                        labels = function(x) sub("^0\\.", ".", as.character(x)),
                        guide = guide_colorbar(order = 2,
                                                  direction = 'horizontal',
                                                  title.position = 'left')) +
    new_scale_fill() +
    geom_tile(data    = dummy %>% filter(y == 3),
              mapping = aes(x = x, y = y, fill = v)) +
    scale_fill_gradient(low      = '#FFFFFF',
                        high     = ccolor,
                        limits   = lim_list[[cctrl]][['dFC']],
                        na.value = NA,
                        name     = 'dFC',
                        labels   = function(x) sub("^0\\.", ".", as.character(x)),
                        guide    = guide_colorbar(order          = 3,
                                                  direction      = 'horizontal',
                                                  title.position = 'left')) +
    theme_void() +
    theme(legend.position = 'top',
          legend.justification = 'center',
          legend.box = 'vertical',
          legend.box.just = 'right',
          legend.text = element_text(size = 8,angle=90),
          legend.title = element_text(size = 9,hjust=1,vjust=0.8,
                                      face='bold',family='Times New Roman'),
          legend.margin = margin(t = 5, b = 5, unit = 'pt'),
          legend.spacing.y = unit(-4, 'pt'),
          plot.margin = grid::unit(c(0, 0, 0, -7), 'mm'))
  plist[[pidx]] <- p_leg
  pidx <- pidx + 1
}

# Reorder from column-first to row-first (state x ctrl).
n_col_entries <- nstates_plot + 1
plist_reordered <- vector('list', nctrl * n_col_entries)
for (ctrlidx in 1:nctrl) {
  for (stidx in 1:n_col_entries) {
    new_idx <- (stidx - 1) * nctrl + ctrlidx
    old_idx <- (ctrlidx - 1) * n_col_entries + stidx
    plist_reordered[[new_idx]] <- plist[[old_idx]]
  }
}

# Save.
outfile <- paste0(outpath, 'allbrain_controlplot.jpg')
ggsave(outfile,
       arrangeGrob(grobs = plist_reordered,
                   nrow = n_col_entries,
                   ncol = nctrl,
                   padding = unit(0,'mm')),
       units = 'in', width = xwidth, height = yheight,
       dpi = 1440, bg = 'white')
