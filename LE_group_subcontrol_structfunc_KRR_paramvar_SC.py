# -*- coding: utf-8 -*-
"""
For a given number of repetitions, inner loop number, and outer loop number for CV, 
find out whether the main findings of the paper hold across different parameter variations
of threshold type and level for the SC model - controllability and strength exhibit high correspondence, and 
g is predicted better than gF and gC. Generate line plots that display these comparisons.
Output:
threshold_C_vs_D.png Compare controllability and strength across thresholds.
threshold_g_vs_gF_vs_gC.png Compare g and gF and gC across thresholds.

Usage: 
    LE_group_subcontrol_structfunc_KRR_paramvar_SC.py <nrep> <inner_k> <outer_k>
    
Arguments:
    
    <nrep> Number of CV repetitions
    <inner_k> Inner K in K-fold CV hyperparameter search
    <outer_k> Outer K in K-fold CV

"""

import os, sys, time, random, colorsys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
from docopt import docopt

if __name__ == '__main__':
    __spec__ = None
    
     #Catches arguments.
    args = docopt(__doc__)
    nrep = args['<nrep>']
    inner_k = args['<inner_k>']
    outer_k = args['<outer_k>']
    print('Doing:',nrep,inner_k,outer_k)

    #Set paths.
    subgroup = 'full'
    sc_subgroup = 'dr_full'
    outpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
             subgroup+'/SC_comkthres/'+sc_subgroup+'/KRRXFS/'+
             nrep+'_'+inner_k+'_'+outer_k+'/')
    os.makedirs(outpath,exist_ok=True)

    # ------------------------ Strength vs Controllability ----------------------- #

    # Threshold compare.
    inkey = ('/rep_testacc')
    sctype = 'fpt'
    k = '6'
    ctrltypes = ['ave','mod','abs_deg'] 
    nctrltypes = len(ctrltypes)
    thresvals = ['10','20','30','40','50','60','70','80','90'] 
    nthresvals = len(thresvals)
    consist_labs = ['none_0'] + [('groupconsist_'+x) for x in thresvals]
    str_labs = ['none_0'] + [('groupstr_'+x) for x in thresvals]
    threslabs = consist_labs + [('groupstr_'+x) for x in thresvals] 
    nthres = len(threslabs)
    statetype = 'SC'
    septype = 'comCFAng'
    cog_septype ='gCFA'
    nrows = nthres
    ncols = nctrltypes

    # Read.
    allmat = pd.DataFrame(np.zeros((nrows, ncols)), index=threslabs, columns=ctrltypes)
    for thidx in range(nthres):
        cthres = threslabs[thidx]
        threstype,thresval = cthres.split('_')      
        for ctrltype in ctrltypes:
            inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                + subgroup + '/'
                + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                + threstype + '/' + thresval + '/KRRXFS/'
                + ctrltype + '_' + statetype + '_' + septype + '_'
                + nrep + '_' + inner_k + '_'
                + outer_k + '_' + sctype)
            infile = (inpath + '/score_collect.h5')
            store = pd.HDFStore(infile, 'r')
            testacc = store.select(inkey)
            store.close()
            allmat.loc[cthres,ctrltype] = testacc.mean(axis=0).loc[cog_septype]
    
    # Plot.
    ax_titles = {
        'groupconsist': 'Consistency',
        'groupstr':     'Strength',
    }
    ctrltype_labels = {
        'ave':     'AC',
        'mod':     'MC',
        'abs_deg': 'S',
    }
    colors = {
        'ave':     '#1f77b4',
        'mod':     '#d62728',
        'abs_deg': '#2ca02c',
    }
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    fig.subplots_adjust(wspace=0.08, right=0.78)
    x_shared = np.arange(len(thresvals) + 1)   
    for ax, group_labs, group_key in [(ax1, consist_labs, 'groupconsist'),
                                      (ax2, str_labs,     'groupstr')]:
        for ctrltype in ctrltypes:
            col = colors[ctrltype]
            y   = np.array([allmat.loc[tl, ctrltype] for tl in group_labs], dtype=float)
            ax.plot(x_shared, y, color=col, linestyle='-', linewidth=2,
                    marker='o', markersize=4, label=ctrltype_labels[ctrltype])
        tick_labels = ['0'] + thresvals
        ax.set_xticks(x_shared)
        ax.set_xticklabels(tick_labels, fontsize=8)
        ax.set_xlabel('Threshold (% Removed)')
        ax.set_title(ax_titles[group_key], fontsize=10)
        ax.set_xlim(-0.3, len(thresvals) + 0.3)
    ax1.set_ylabel('R²')
    ctrl_handles = [
        mlines.Line2D([], [], color=colors[ct], linewidth=2,
                      label=ctrltype_labels[ct])
        for ct in ctrltypes
    ]
    fig.legend(handles=ctrl_handles, title='Metric',
               loc='upper left',
               bbox_to_anchor=(0.80, 0.88),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    plt.savefig((outpath+'/threshold_C_vs_D.png'),dpi=720,bbox_inches='tight')
    plt.close()

    # ------------------------------- g vs gF vs gC ------------------------------ #

    # Threshold compare.
    inkey = ('/rep_testacc')
    sctype = 'fpt'
    k = '6'
    ctrltypes = ['ave','mod','abs_deg'] 
    nctrltypes = len(ctrltypes)
    thresvals = ['10','20','30','40','50','60','70','80','90'] 
    nthresvals = len(thresvals)
    consist_labs = ['none_0'] + [('groupconsist_'+x) for x in thresvals]
    str_labs = ['none_0'] + [('groupstr_'+x) for x in thresvals]
    threslabs = consist_labs + [('groupstr_'+x) for x in thresvals] 
    nthres = len(threslabs)
    statetype = 'SC'
    septype = 'comCFAng'
    cog_septype = ['gCFA','P24_CR','PV'] 
    ncog = len(cog_septype)
    nrows = nthres
    ncols = ncog

    # Read.
    avemat = pd.DataFrame(np.zeros((nrows,ncols)),index=threslabs,columns=cog_septype)
    modmat = pd.DataFrame(np.zeros((nrows,ncols)),index=threslabs,columns=cog_septype)
    degmat = pd.DataFrame(np.zeros((nrows,ncols)),index=threslabs,columns=cog_septype)
    for thidx in range(nthres):
        cthres = threslabs[thidx]
        threstype,thresval = cthres.split('_')
        for cogidx in range(ncog):
            ccog = cog_septype[cogidx]
            for ctrltype, outmat in [('ave', avemat), ('mod', modmat), ('abs_deg', degmat)]:
                inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                    + subgroup + '/'
                    + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                    + threstype + '/' + thresval + '/KRRXFS/'
                    + ctrltype + '_' + statetype + '_' + septype + '_'
                    + nrep + '_' + inner_k + '_'
                    + outer_k + '_' + sctype)
                infile = (inpath + '/score_collect.h5')
                store = pd.HDFStore(infile, 'r')
                testacc = store.select(inkey)
                store.close()
                outmat.loc[cthres,ccog] = testacc.mean(axis=0).loc[ccog]
    
    # Plot.
    ax_titles = {
        'groupconsist': 'Consistency',
        'groupstr':     'Strength',
    }
    cog_labels = {
        'gCFA':    'g',
        'P24_CR':  'gF',
        'PV':      'gC',
    }
    ctrltype_labels = {
        'ave': 'AC',
        'mod': 'MC',
        'abs_deg': 'S'
    }
    cog_colors = {
        'gCFA':    '#1f77b4',
        'P24_CR':  '#d62728',
        'PV':      '#2ca02c',
    }
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    fig.subplots_adjust(wspace=0.08, right=0.78)
    x_shared = np.arange(len(thresvals) + 1)
    for ax, group_labs, group_key in [(ax1, consist_labs, 'groupconsist'),
                                      (ax2, str_labs,     'groupstr')]:
        for ccog in cog_septype:
            col   = cog_colors[ccog]
            ave_y = np.array([avemat.loc[tl, ccog] for tl in group_labs], dtype=float)
            mod_y = np.array([modmat.loc[tl, ccog] for tl in group_labs], dtype=float)
            deg_y = np.array([degmat.loc[tl, ccog] for tl in group_labs], dtype=float)
            ax.plot(x_shared, ave_y, color=col, linestyle='-',  linewidth=2, marker='o', markersize=4)
            ax.plot(x_shared, mod_y, color=col, linestyle='--', linewidth=2, marker='o', markersize=4)
            ax.plot(x_shared, deg_y, color=col, linestyle=':', linewidth=2, marker='o', markersize=4)
        tick_labels = ['0'] + thresvals
        ax.set_xticks(x_shared)
        ax.set_xticklabels(tick_labels, fontsize=8)
        ax.set_xlabel('Threshold (% Removed)')
        ax.set_title(ax_titles[group_key], fontsize=10)
        ax.set_xlim(-0.3, len(thresvals) + 0.3)
    ax1.set_ylabel('R²')
    cog_handles = [
        mlines.Line2D([], [], color=cog_colors[cg], linewidth=2,
                      label=cog_labels[cg])
        for cg in cog_septype
    ]
    style_handles = [
        mlines.Line2D([], [], color='grey', linestyle='-',  linewidth=2,
                      label=ctrltype_labels['ave']),
        mlines.Line2D([], [], color='grey', linestyle='--', linewidth=2,
                      label=ctrltype_labels['mod']),
        mlines.Line2D([], [], color='grey', linestyle=':', linewidth=2,
                      label=ctrltype_labels['abs_deg'])              
    ]
    fig.legend(handles=cog_handles, title='Domain',
               loc='upper left',
               bbox_to_anchor=(0.80, 0.92),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    fig.legend(handles=style_handles, title='Metric',
               loc='upper left',
               bbox_to_anchor=(0.80, 0.35),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    plt.savefig((outpath+'/threshold_g_vs_gF_vs_gC.png'), dpi=720, bbox_inches='tight')
    plt.close()
    print('Parameter variation done.')
