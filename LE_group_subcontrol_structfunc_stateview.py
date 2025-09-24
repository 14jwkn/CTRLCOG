#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the given k for clustering, SC type of normalization, and thresholding type
and percentage, plot the dFC, SC, and sFC matrices organized according to
hemisphere and network. Also save the network legends and color bars.
Output:
dFC_netleg.jpg dFC states network legend.
dFC_colleg.jpg dFC states color bar. 
'dFC_S'+str(kidx+1)+'_scaled.jpg' dFC state matrix plotted.
sctype+'_SC_unscaled.jpg' SC matrix plotted.
sctype+'_SC_netleg.jpg' SC network legend.
sctype+'_SC_colleg.jpg' SC color bar.
sFC_unscaled.jpg sFC matrix plotted.
sFC_netleg.jpg sFC network legend.
sFC_colleg.jpg sFC color bar.

Usage: 
    LE_group_subcontrol_structfunc_stateview.py <k> <sctype> <threstype> <thresval>
    
Arguments:
    
    <k> k for k-means clustering
    <sctype> SC type of normalization
    <threstype> Threshold type
    <thresval> Threshold value

"""

import os, h5py
import numpy as np
import pandas as pd
import seaborn as sns
import scicomap as sc
from docopt import docopt
from collections import OrderedDict
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable

#Convert vectorized into sorted matrix.
def visvec(dFCwin,netlabels):
    
    #Format.
    corrmat = pd.DataFrame(dFCwin)
    
    #Add a column for the labels.
    rowlabelled = pd.concat([pd.Series(netlabels),pd.DataFrame(corrmat)],axis=1)
    
    #Add a row for the labels.
    colnetlabels = [0] + netlabels
    rowlabelled.loc[-1] = colnetlabels
    rowlabelled.index = rowlabelled.index + 1
    collabelled = rowlabelled.sort_index()
    collabelled.columns = range(361)
    
    #Adds axes labels to enable reference.
    collabelled = collabelled.rename_axis('Index')
    collabelled = collabelled.rename_axis('Columns',axis='columns')
    
    #Sort the rows and columns.
    rowsort = collabelled.sort_values(by=[0,'Index'],axis=0)
    colsort = rowsort.sort_values(by=[0,'Columns'],axis=1)
    
    #Reset indices. Save the matrix to list.
    reformatted = colsort.reset_index(drop=True)
    reformatted.columns = range(reformatted.shape[1])
    return reformatted

#Convert vectorized into sorted matrix by hemisphere.
def visvec_hemi(dFCwin,netlabels):
    
    #Format.
    corrmat = pd.DataFrame(dFCwin)
    
    #Add a column for the network labels.
    rowlabelled = pd.concat([pd.Series(netlabels),pd.DataFrame(corrmat)],axis=1)
    
    #Add a row for the network labels.
    colnetlabels = [0] + netlabels
    rowlabelled.loc[-1] = colnetlabels
    rowlabelled.index = rowlabelled.index + 1
    collabelled = rowlabelled.sort_index()
    collabelled.columns = range(361)

    #Add a column for hemisphere labels.
    hemilabels = [0] + [1]*180 + [2]*180
    rowlabelled = pd.concat([pd.Series(hemilabels),collabelled],axis=1)

    #Add a row for the network labels.
    colhemilabels = [0] + hemilabels
    rowlabelled.loc[-1] = colhemilabels
    rowlabelled.index = rowlabelled.index + 1
    collabelled = rowlabelled.sort_index()
    collabelled.columns = range(362)
    
    #Adds axes labels to enable reference.
    collabelled = collabelled.rename_axis('Index')
    collabelled = collabelled.rename_axis('Columns',axis='columns')
    
    #Sort the rows and columns.
    rowsort = collabelled.sort_values(by=[0,1,'Index'],axis=0)
    colsort = rowsort.sort_values(by=[0,1,'Columns'],axis=1)
    
    #Reset indices. Save the matrix to list.
    reformatted = colsort.reset_index(drop=True)
    reformatted.columns = range(reformatted.shape[1])
    return reformatted

if __name__ == '__main__':
    __spec__ = None
    
    #Catches arguments.
    args = docopt(__doc__)
    k = args['<k>']
    sctype = args['<sctype>']
    threstype = args['<threstype>']
    thresval = args['<thresval>']
    print('Doing:',k,sctype,threstype,thresval)
    
    #Set up I/O.
    subgroup = 'full'
    sc_subgroup = 'dr_full'
    basepath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                subgroup+'/'+k+'/SC_dFC/'+sc_subgroup+'/collect/'+
                threstype+'/'+thresval+'/state_images/')
    
    #Set general parameters.
    nk = int(k)
    
    #Set parameters for matrix plots.
    labsize = 18
    negcode = 250
    poscode = 20
    satcode = 100
    lumcode = 60
  
    #Read in network labels and convert them to integers.
    with open('colenetlabels.txt') as f:
        netlabels = [label.rstrip() for label in f] 
    netlabels = list(map(int,netlabels))
    
    #Generate colorbar components from the network labels.
    groups = sorted(netlabels)
    groupstr = []
    namelist = ['Primary Visual (VIS1)','Secondary Visual (VIS2)','Somatomotor (SMN)',
                'Cingulo-Opercular (CON)','Dorsal Attention (DAN)','Language (LAN)',
                'Frontoparietal (FPN)','Auditory (AUD)','Default Mode (DMN)',
                'Posterior Multimodal (PMM)','Ventral Multimodal (VMM)',
                'Orbitoaffective (ORA)']
    for net in groups:
        groupstr.append(namelist[int(net)-1])
    lut_colors = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b',
                  '#e377c2','#7f7f7f','#bcbd22','#17becf','#AEFD6C','#fff200']
    lut_dict = OrderedDict(zip(namelist,lut_colors))
    cmap = colors.LinearSegmentedColormap.from_list('Networks',lut_colors,N=len(namelist))
    class_labels = np.array(groups).reshape((1,360))
    class_labels = np.subtract(class_labels,1)

    #Read in the MNI coordinates.
    glasscoords = np.loadtxt('glasser_coords.txt') 

    #Generate color names for network labels.
    net_labels = []
    for net in netlabels:
        netstr = namelist[int(net)-1]
        net_labels.append(lut_dict[netstr])
    net_labels = np.array(net_labels)

    #Read in the dFC states and limits.
    inpath = basepath
    dFC_states = []
    dFC_limvals = []
    for i in range(int(k)):

        #Read in the states.
        infile = (inpath+'SC_sFC_dFC.h5')
        instore = h5py.File(infile,'r')
        inkey = ('/s'+str(i+1))
        inmat = np.array(instore[inkey]).T
        instore.close()
        dFC_states.append(inmat)
        inmat = np.copy(inmat)
        inmat[np.diag_indices_from(inmat)] = float('nan')
        dFC_limvals.append(np.nanmin(inmat))
        dFC_limvals.append(np.nanmax(inmat))  
      
    #Reformat each state.
    dFC_collect = []
    for i in range(int(k)):
        
        #Select the state, reformat, and append.
        print('Doing state:',i+1)
        dFCwin = dFC_states[i]
        reformatted = visvec(dFCwin,netlabels)
        dFC_collect.append(reformatted.iloc[1:,1:].values)
        new_class_labels = reformatted.iloc[0,1:].values
        new_class_labels = np.expand_dims(new_class_labels,axis=0)
    
    #Find limits.
    dFCmin = min(dFC_limvals)   
    dFCmax = max(dFC_limvals)       
    
    #Save legend.
    handles = [Patch(facecolor=lut_dict[name]) for name in lut_dict]
    plt.legend(handles,lut_dict,title='Networks',loc='center') 
    plt.axis('off')
    outfile = basepath+'dFC_netleg.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close()  
    
    #Save color bar.
    if dFCmin != 0:
        divnorm=colors.TwoSlopeNorm(vmin=dFCmin,vcenter=0,vmax=dFCmax)
    else:
        divnorm=colors.TwoSlopeNorm(vcenter=0,vmax=dFCmax)
    c = plt.pcolormesh(dFC_collect[0],cmap='RdBu_r',norm=divnorm)
    cbar = plt.colorbar(c,spacing='proportional')
    cbar.ax.tick_params(labelsize=labsize)
    cbar.ax.set_yscale('linear')
    outfile = basepath+'dFC_colleg.jpg'
    plt.savefig(outfile,dpi=720,bbox_inches='tight')
    plt.close()  

    #Plot each dFC state individually for higher resolution.
    for kidx in range(nk):
        
        #Get state data.
        fig, ax = plt.subplots()
        rawdat = dFC_collect[kidx]

        #Make the diagonal the max value.
        dat = rawdat.copy()
        dat[np.diag_indices_from(dat)] = np.nanmax(dat)

        #Plot matrix.
        if dFCmin != 0:
            divnorm=colors.TwoSlopeNorm(vmin=dFCmin,vcenter=0,vmax=dFCmax)
        else:
            divnorm=colors.TwoSlopeNorm(vcenter=0,vmax=dFCmax)
        im = ax.imshow(dat,cmap='RdBu_r',norm=divnorm)
        
        #Plot network bars.
        divider = make_axes_locatable(ax)
        x_ax = divider.append_axes('top',size='10%',pad=0)
        y_ax = divider.append_axes('left',size='10%',pad=0)
        x_ax.imshow(new_class_labels,aspect='auto',cmap=cmap)
        y_ax.imshow(np.transpose(new_class_labels),aspect='auto',cmap=cmap)
        
        #Remove axes values.
        ax.axis('off')
        x_ax.set_axis_off()
        y_ax.set_axis_off()    

        #Save.
        outfile = basepath+'dFC_S'+str(kidx+1)+'_scaled.jpg'
        plt.savefig(outfile,dpi=720)
        plt.close()
        
    #Read in the SC state and produced sorted version.
    inpath = basepath
    infile = (inpath+'SC_sFC_dFC.h5')
    instore = h5py.File(infile,'r')
    inkey = ('/SC_'+sctype)
    scmat = np.array(instore[inkey]).T
    instore.close()
    inmat = np.copy(scmat)
    inmat[np.diag_indices_from(inmat)] = float('nan') 
    reformatted = visvec_hemi(scmat,netlabels)
    new_class_labels = reformatted.iloc[1,2:].values
    new_class_labels = np.subtract(new_class_labels,1)
    new_class_labels = np.expand_dims(new_class_labels,axis=0)
      
    #Convert zeros to a low number, take the log, and add
    #by a constant to make the minimum value zero. Repeat for the 
    #unsorted version.
    print('Doing SC')
    sc_view = reformatted.iloc[2:,2:].values
    sc_view[np.where(sc_view==0)]=0.000001
    sc_view = np.log(sc_view)
    sc_view = np.add(sc_view,-np.nanmin(sc_view))
    sc_min = np.nanmin(sc_view)
    sc_max = np.nanmax(sc_view) 
    raw_sc_view = inmat.copy()
    raw_sc_view[np.where(raw_sc_view==0)]=0.000001
    raw_sc_view = np.log(raw_sc_view)
    raw_sc_view = np.add(raw_sc_view,-np.nanmin(raw_sc_view))
    
    #Plot for SC, unscaled.
    fig, ax = plt.subplots()
    rawdat = sc_view

    #Make the diagonal the max value.
    dat = rawdat.copy()
    dat[np.diag_indices_from(dat)] = np.nanmax(dat)

    #Plot matrix.
    divnorm = colors.TwoSlopeNorm(vcenter=0)
    im = ax.imshow(dat,cmap='RdBu_r',norm=divnorm)

    #Plot network bars.
    divider = make_axes_locatable(ax)
    x_ax = divider.append_axes('top',size='10%',pad=0)
    y_ax = divider.append_axes('left',size='10%',pad=0)
    x_ax.imshow(new_class_labels,aspect='auto',cmap=cmap)
    y_ax.imshow(np.transpose(new_class_labels),aspect='auto',cmap=cmap)
        
    #Remove axes values.
    ax.axis('off')
    x_ax.set_axis_off()
    y_ax.set_axis_off()    
    
    #Save.
    outfile = basepath+sctype+'_SC_unscaled.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close()      
    
    #Save legend.
    handles = [Patch(facecolor=lut_dict[name]) for name in lut_dict]
    plt.legend(handles,lut_dict,title='Networks',loc='center') 
    plt.axis('off')
    outfile = basepath+sctype+'_SC_netleg.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close()  
    
    #Save color bar.
    divnorm=colors.TwoSlopeNorm(vcenter=0)
    c = plt.pcolormesh(sc_view,cmap='RdBu_r',norm=divnorm)
    cbar = plt.colorbar(c)
    cbar.ax.tick_params(labelsize=labsize)
    outfile = basepath+sctype+'_SC_colleg.jpg'
    plt.savefig(outfile,dpi=720,bbox_inches='tight')
    plt.close()  
    
    #Read in the sFC state and limits.
    inpath = basepath
    infile = (inpath+'SC_sFC_dFC.h5')
    instore = h5py.File(infile,'r')
    inkey = ('/sFC')
    sFCmat = np.array(instore[inkey]).T
    instore.close()
    inmat = np.copy(sFCmat)
    inmat[np.diag_indices_from(inmat)] = float('nan') 
    sFC_min = np.nanmin(sFCmat)
    sFC_max = np.nanmax(sFCmat) 
      
    #Reformat.
    print('Doing sFC')
    reformatted = visvec(sFCmat,netlabels)
    sFC_view = reformatted.iloc[1:,1:].values
    new_class_labels = reformatted.iloc[0,1:].values
    new_class_labels = np.expand_dims(new_class_labels,axis=0)

    #Plot for sFC, unscaled.
    fig, ax = plt.subplots()
    rawdat = sFC_view

    #Make the diagonal the max value.
    dat = rawdat.copy()
    dat[np.diag_indices_from(dat)] = np.nanmax(dat)

    #Plot matrix.
    divnorm = colors.TwoSlopeNorm(vcenter=0)
    im = ax.imshow(dat,cmap='RdBu_r',norm=divnorm)

    #Plot network bars.
    divider = make_axes_locatable(ax)
    x_ax = divider.append_axes('top',size='10%',pad=0)
    y_ax = divider.append_axes('left',size='10%',pad=0)
    x_ax.imshow(new_class_labels,aspect='auto',cmap=cmap)
    y_ax.imshow(np.transpose(new_class_labels),aspect='auto',cmap=cmap)
        
    #Remove axes values.
    ax.axis('off')
    x_ax.set_axis_off()
    y_ax.set_axis_off()    
    
    #Save.
    outfile = basepath+'sFC_unscaled.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close()      
    
    #Save legend.
    handles = [Patch(facecolor=lut_dict[name]) for name in lut_dict]
    plt.legend(handles,lut_dict,title='Networks',loc='center') 
    plt.axis('off')
    outfile = basepath+'sFC_netleg.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close()  

    #Save color bar.
    divnorm=colors.TwoSlopeNorm(vcenter=0)
    c = plt.pcolormesh(sFC_view,cmap='RdBu_r',norm=divnorm)
    cbar = plt.colorbar(c,spacing='proportional')
    cbar.ax.tick_params(labelsize=labsize)
    cbar.ax.set_yscale('linear')
    outfile = basepath+'sFC_colleg.jpg'
    plt.savefig(outfile,dpi=720,bbox_inches='tight')
    plt.close() 
