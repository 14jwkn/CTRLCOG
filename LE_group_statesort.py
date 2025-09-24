#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the set of LE or dFC centroids across a range of k, print out the 2D images and 
generate a table which sorts them. Then, visually sort them and edit the table. 
Run this code again with the edited table to visualize the sorted plots.
Output:
subgroup+'_sortclass.csv' Sorting matrix containing the original ID and new ID.
'maxmin_k'+statetype+'.jpg' 2D plots where the maximum and minimum values are cross-k.
'legend_k'+statetype+'.jpg' Legend for the network labels.
'maxmin_clbr_k'+statetype+'_noticks.jpg' Color bar for the max-min 2D plots with no ticks.
'maxmin_clbr_k'+statetype+'_withticks.jpg' Color bar for the max-min 2D plots with ticks.
'default_k'+statetype+'.jpg' 2D plots with automatic max-min for each plot.
'k'+statetype+'.jpg' 2D plots with no labelling.

Usage: 
    LE_group_statesort.py <statetype> 
    
Arguments:

    <statetype> dFC or LE

"""

import os, h5py
import numpy as np
import pandas as pd
import seaborn as sns
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
    
    #Sort the rows and columns.w
    rowsort = collabelled.sort_values(by=[0,'Index'],axis=0)
    colsort = rowsort.sort_values(by=[0,'Columns'],axis=1)
    
    #Reset indices. Save the matrix to list.
    reformatted = colsort.reset_index(drop=True)
    reformatted.columns = range(reformatted.shape[1])
    return reformatted

if __name__ == '__main__':
    __spec__ = None
    
    #Catches command-line arguments.
    args = docopt(__doc__)
    statetype = args['<statetype>']
    print('Doing:',statetype)
    
    #Set up I/O.
    subfile = 'r_full_submain.txt'
    subgroup = 'full'
    outpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/visclass/'+
              subgroup+'/allk/')                   
    os.makedirs(outpath,exist_ok=True)
        
    #Read in network labels and convert them to integers.
    with open('colenetlabels.txt') as f:
        netlabels = [label.rstrip() for label in f] 
    netlabels = list(map(int,netlabels))
    
    #Generate colorbar components from the network labels.
    groups = sorted(netlabels)
    groupstr = []
    namelist = ['Primary Visual (VIS1)','Secondary Visual (VIS2)',
                'Somatomotor (SMN)','Cingulo-Opercular (CON)',
                'Dorsal Attention (DAN)','Language (LAN)',
                'Frontoparietal (FPN)','Auditory (AUD)','Default Mode (DMN)',
                'Posterior Multimodal (PMM)','Ventral Multimodal (VMM)',
                'Orbitoaffective (ORA)']
    for net in groups:
        groupstr.append(namelist[int(net)-1])
    lut_colors = list(sns.color_palette(palette='tab10',n_colors=len(namelist)))
    lut_dict = OrderedDict(zip(namelist,sns.color_palette(palette='tab10',n_colors=len(namelist))))
    cmap = colors.LinearSegmentedColormap.from_list('Networks',lut_colors,N=len(namelist))
    class_labels = np.array(groups).reshape((1,360))
    class_labels = np.subtract(class_labels,1)
    
    #Set desired ks.
    kstart = 2
    kend = 12
    knum = kend - kstart + 1
    klist = [str(x) for x in range(kstart,kend+1)]
    
    #Read file for rearrangement or create it if it does not exist.
    infile = (outpath+subgroup+'_sortclass.csv')
    if os.path.exists(infile):
        reorder_mat = pd.read_csv(infile,header=0,index_col=0)
        reorder_mat.index = [str(x) for x in reorder_mat.index]
        reorder_mat.columns = [str(x) for x in reorder_mat.columns] 
    else:   
        reorder_mat = pd.DataFrame(np.full((knum,kend),''))
        reorder_mat.index = klist
        reorder_mat.columns = [str(x) for x in range(1,kend+1)]
        reorder_mat.to_csv(infile,index=True,header=True)
        reorder_mat = reorder_mat.replace('',None)
    
    #For each k.
    state_collect = []
    limvals = []
    for kidx in range(knum):
        
        #Extract.
        k = klist[kidx]
        
        #Read in the states.
        inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                  subgroup+'/'+k+'/') 
        states_dFC = []
        for i in range(int(k)):
            
            #Read in the states.
            infile = (inpath+'k_'+statetype+'.h5')
            instore = h5py.File(infile,'r')
            inkey = ('/'+subgroup+'_'+str(i+1))
            inmat = np.array(instore[inkey]).T
            states_dFC.append(inmat)
            
            #Find limits.
            inmat = np.copy(inmat)
            inmat[np.diag_indices_from(inmat)] = float('nan')
            limvals.append(np.nanmin(inmat))
            limvals.append(np.nanmax(inmat))
            instore.close()
        
        #Read the shuffler and shuffle if it the order has been picked.
        reorder_k = reorder_mat.loc[k,:].dropna().values.astype(int).tolist()
        if reorder_k != []:
            states_dFC = [states_dFC[i-1] for i in reorder_k]
      
        #For each state.
        for i in range(int(k)):
            
            #Select the state, reformat, and append.
            print('Doing state:',i+1)
            dFCwin = states_dFC[i]
            reformatted = visvec(dFCwin,netlabels)
            state_collect.append(reformatted.iloc[1:,1:].values)
        
        #Add empty parts in remaining ones.
        kdiff = kend - int(k)
        for i in range(kdiff):
            state_collect.append(np.zeros((360,360)))
   
    #Find limits.
    dFCmin = min(limvals)   
    dFCmax = max(limvals)     
    
    #Generate labels.
    plabs = []
    for kidx in range(knum):
        k = klist[kidx]
        for i in range(int(k)):
            plabs.append(str(i+1))        
        
    #Plot each k with maxmin values.
    pidx = 0
    fig, axs = plt.subplots(nrows=knum,ncols=kend)
    for ax, dat in zip(axs.ravel(),state_collect):
        if np.all(dat == 0):
            ax.set_visible(False)
        else:
            
            #Create color bar axes.
            divider = make_axes_locatable(ax)
            x_ax = divider.append_axes('top',size='10%',pad=0)
            y_ax = divider.append_axes('left',size='10%',pad=0)
            divnorm=colors.TwoSlopeNorm(vmin=dFCmin,vcenter=0,vmax=dFCmax)
            
            #Plot.
            im = ax.imshow(dat,cmap='RdBu_r',norm=divnorm)
            x_ax.imshow(class_labels,aspect='auto',cmap=cmap)
            y_ax.imshow(np.transpose(class_labels),aspect='auto',cmap=cmap)
            
            #Remove axes values.
            ax.axis('off')
            x_ax.set_axis_off()
            y_ax.set_axis_off()    
            
            #Add title.
            pstate = plabs[pidx]
            ax.set_title(pstate,fontweight='bold',loc='left',x=-0.25,y=0.77,size=3.5)
            
            #Increment.
            pidx = pidx + 1
            
    #Save.
    plt.subplots_adjust(wspace=0.2,hspace=0.09,
                        top=0.99,bottom=0.01,left=0.10,right=0.9)
    outfile = outpath+'maxmin_k'+statetype+'.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close()    
    
    #Save legend.
    handles = [Patch(facecolor=lut_dict[name]) for name in lut_dict]
    plt.legend(handles,lut_dict,title='Networks',loc='center') 
    plt.axis('off')
    outfile = outpath+'legend_k'+statetype+'.jpg'
    plt.savefig(outfile,dpi=720,bbox_inches='tight')
    plt.close()  
    
    #Save color bar.
    plt.imshow(state_collect[0],cmap='RdBu_r',vmin=dFCmin,vmax=dFCmax)    
    cbar = plt.colorbar()
    cbar.ax.tick_params(size=0,labelsize=15) 
    outfile = outpath+'maxmin_clbr_k'+statetype+'_noticks.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close() 
    plt.imshow(state_collect[0],cmap='RdBu_r',vmin=dFCmin,vmax=dFCmax)    
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=15) 
    outfile = outpath+'maxmin_clbr_k'+statetype+'_withticks.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close() 
    
    #Plot each k without maxmin values.
    pidx = 0
    fig, axs = plt.subplots(nrows=knum,ncols=kend)
    for ax, dat in zip(axs.ravel(),state_collect):
        if np.all(dat == 0):
            ax.set_visible(False)
        else:
            
            #Create color bar axes.
            divider = make_axes_locatable(ax)
            x_ax = divider.append_axes('top',size='10%',pad=0)
            y_ax = divider.append_axes('left',size='10%',pad=0)
            divnorm=colors.TwoSlopeNorm(vcenter=0)
            
            #Plot.
            im = ax.imshow(dat,cmap='RdBu_r',norm=divnorm)
            x_ax.imshow(class_labels,aspect='auto',cmap=cmap)
            y_ax.imshow(np.transpose(class_labels),aspect='auto',cmap=cmap)
            
            #Remove axes values.
            ax.axis('off')
            x_ax.set_axis_off()
            y_ax.set_axis_off()    
            
            #Add title.
            pstate = plabs[pidx]
            ax.set_title(pstate,fontweight='bold',loc='left',x=-0.25,y=0.77,size=3.5)
            
            #Increment.
            pidx = pidx + 1
            
    #Save.
    plt.subplots_adjust(wspace=0.2,hspace=0.09,
                        top=0.99,bottom=0.01,left=0.10,right=0.9)
    outfile = outpath+'default_k'+statetype+'.jpg'
    plt.savefig(outfile,dpi=720)
    plt.close()    
    
    #Plot each k with no labelling.
    fig, axs = plt.subplots(nrows=knum,ncols=kend)
    for ax, dat in zip(axs.ravel(),state_collect):
        if np.all(dat == 0):
            ax.set_visible(False)
        else:
            
            #Plot.
            im = ax.imshow(dat,cmap='RdBu_r',norm=divnorm)
            
            #Remove axes values.
            ax.axis('off')
              
    #Save.
    outfile = outpath+'k'+statetype+'.jpg'
    plt.savefig(outfile,bbox_inches='tight',dpi=720)
    plt.close() 
            
