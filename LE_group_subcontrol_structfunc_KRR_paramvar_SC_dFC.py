# -*- coding: utf-8 -*-
"""
For a given number of repetitions, inner loop number, and outer loop number for CV, 
find out whether the main findings of the paper hold across different parameter variations
for the SC & dFC model - dFC performs better than SC and controllability and strength exhibit
high correspondence, g is predicted better than gF and gC, the principal gradient exhibits a relationship 
with regional importance in the expected direction for dFC. The parameters that are varied include threshold 
type and level, k for dFC, cognitive variable measures, and SC normalization. Generate line plots 
that display these comparisons.
Output:
threshold_SC_vs_dFC.png Compare dFC and SC, and controllability and strength, across thresholds.
k_SC_vs_dFC.png Compare dFC and SC, and controllability and strength, across k.
altcog_SC_vs_dFC.png Compare dFC and SC, and controllability and strength, across cognitive variable measures.
scnorm_SC_vs_dFC.png Compare dFC and SC, and controllability and strength, across SC normalization types.
threshold_g_vs_gF_vs_gC.png Compare g and gF and gC across thresholds.
k_g_vs_gF_vs_gC.png Compare g and gF and gC across k.
altcog_g_vs_gF_vs_gC.png Compare g and gF and gC across cognitive variable measures.
scnorm_g_vs_gF_vs_gC.png Compare g and gF and gC across SC normalization types.
threshold_sFCgr.png Compare the principal gradient and regional importance across thresholds.
k_sFCgr.png Compare the principal gradient and regional importance across k.
altcog_sFCgr.png Compare the principal gradient and regional importance across cognitive variable measures.
scnorm_sFCgr.png Compare the principal gradient and regional importance across SC normalization types.

Usage: 
    LE_group_subcontrol_structfunc_KRR_paramvar_SC_dFC.py <nrep> <inner_k> <outer_k>
    
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
             subgroup+'/SC_dFC_comkthres/'+sc_subgroup+'/KRRXFS/'+
             nrep+'_'+inner_k+'_'+outer_k+'/')
    os.makedirs(outpath,exist_ok=True)

    # --------------------------------- SC vs dFC -------------------------------- #

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
    statetypes = ['SC','dFCcat','SC_dFCcat'] 
    nstatetypes = len(statetypes)
    septype = 'comCFAng'
    cog_septype ='gCFA'
    nrows = nthres
    ncols = nstatetypes

    # Read.
    avemat = pd.DataFrame(np.zeros((nrows,ncols)),index=threslabs,columns=statetypes)
    modmat = pd.DataFrame(np.zeros((nrows,ncols)),index=threslabs,columns=statetypes)
    degmat = pd.DataFrame(np.zeros((nrows,ncols)),index=threslabs,columns=statetypes)
    for thidx in range(nthres):
        cthres = threslabs[thidx]
        threstype,thresval = cthres.split('_')
        for stidx in range(nstatetypes):
            statetype = statetypes[stidx]
            for ctrltype, outmat in [('ave', avemat), ('mod', modmat), ('abs_deg',degmat)]:
                inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                    + subgroup + '/'
                    + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                    + threstype + '/' + thresval + '/KRRXFS/'
                    + ctrltype + '_' + statetype + '_' + septype + '_'
                    + '_' + nrep + '_' + inner_k + '_'
                    + outer_k + '_' + sctype)
                infile = (inpath + '/score_collect.h5')
                store = pd.HDFStore(infile, 'r')
                testacc = store.select(inkey)
                store.close()
                outmat.loc[cthres,statetype] = testacc.mean(axis=0).loc[cog_septype]
    
    # Plot.
    ax_titles = {
        'groupconsist': 'Consistency',
        'groupstr': 'Strength',
    }
    statetype_labels = {
        'SC':        'SC',
        'dFCcat':    'dFC',
        'SC_dFCcat': 'SC & dFC',
    }
    ctrltype_labels = {
        'ave': 'AC',
        'mod': 'MC',
        'abs_deg': 'S'
    }
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    fig.subplots_adjust(wspace=0.08)
    colors = {
        'SC':        '#1f77b4',
        'dFCcat':    '#d62728',
        'SC_dFCcat': '#2ca02c',
    }
    x_shared = np.arange(len(thresvals) + 1)
    for ax, group_labs, group_key in [(ax1, consist_labs,'groupconsist'),
                                  (ax2, str_labs,'groupstr')]:
        for statetype in statetypes:
            col   = colors[statetype]
            ave_y = np.array([avemat.loc[tl, statetype] for tl in group_labs], dtype=float)
            mod_y = np.array([modmat.loc[tl, statetype] for tl in group_labs], dtype=float)
            deg_y = np.array([degmat.loc[tl, statetype] for tl in group_labs], dtype=float)
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
    color_handles = [
        mlines.Line2D([], [], color=colors[st], linewidth=2,
                    label=statetype_labels[st])
        for st in statetypes
    ]
    style_handles = [
        mlines.Line2D([], [], color='grey', linestyle='-',  linewidth=2,
                    label=ctrltype_labels['ave']),
        mlines.Line2D([], [], color='grey', linestyle='--', linewidth=2,
                    label=ctrltype_labels['mod']),
        mlines.Line2D([], [], color='grey', linestyle=':', linewidth=2,
                    label=ctrltype_labels['abs_deg'])            
    ]
    fig.legend(handles=color_handles, title='Model',
            loc='upper left',
            bbox_to_anchor=(0.92, 0.92), 
            bbox_transform=fig.transFigure,
            borderaxespad=0, framealpha=0.8)
    fig.legend(handles=style_handles, title='Metric',
            loc='upper left',
            bbox_to_anchor=(0.92, 0.45),
            bbox_transform=fig.transFigure,
            borderaxespad=0, framealpha=0.8)
    plt.savefig((outpath+'/threshold_SC_vs_dFC.png'), dpi=720, bbox_inches='tight')
    plt.close()

    # k compare.
    inkey = ('/rep_testacc')
    sctype = 'fpt'
    ks = [str(x) for x in range(2,13)]
    nks = len(ks)
    ctrltypes = ['ave','mod','abs_deg'] 
    nctrltypes = len(ctrltypes)
    thresval = '50'
    threstype = 'groupconsist'
    statetypes = ['SC','dFCcat','SC_dFCcat'] 
    nstatetypes = len(statetypes)
    septype = 'comCFAng'
    cog_septype ='gCFA'
    nrows = nks
    ncols = nstatetypes

    # Read.
    avemat = pd.DataFrame(np.zeros((nrows,ncols)),index=ks,columns=statetypes)
    modmat = pd.DataFrame(np.zeros((nrows,ncols)),index=ks,columns=statetypes)
    degmat = pd.DataFrame(np.zeros((nrows,ncols)),index=ks,columns=statetypes)
    for kidx in range(nks):
        k = ks[kidx]
        for stidx in range(nstatetypes):
            statetype = statetypes[stidx]
            for ctrltype, outmat in [('ave', avemat), ('mod', modmat), ('abs_deg', degmat)]:
                if statetype=='SC':
                    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                            + '/' + subgroup + '/'
                            '6/SC_dFC/' + sc_subgroup + '/collect/'
                            + threstype + '/' + thresval + '/KRRXFS/'
                            + ctrltype + '_' + statetype + '_' + septype + '_'
                            + '_' + nrep + '_' + inner_k + '_'
                            + outer_k + '_' + sctype)
                else:
                    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                        + subgroup + '/'
                        + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                        + threstype + '/' + thresval + '/KRRXFS/'
                        + ctrltype + '_' + statetype + '_' + septype + '_'
                        + '_' + nrep + '_' + inner_k + '_'
                        + outer_k + '_' + sctype)
                infile = (inpath + '/score_collect.h5')
                store = pd.HDFStore(infile, 'r')
                testacc = store.select(inkey)
                store.close()
                outmat.loc[k,statetype] = testacc.mean(axis=0).loc[cog_septype]
    
    # Plot.
    statetype_labels = {
        'SC':        'SC',
        'dFCcat':    'dFC',
        'SC_dFCcat': 'SC & dFC',
    }
    ctrltype_labels = {
        'ave': 'AC',
        'mod': 'MC',
        'abs_deg': 'S'
    }
    colors = {
        'SC':        '#1f77b4',
        'dFCcat':    '#d62728',
        'SC_dFCcat': '#2ca02c',
    }
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.subplots_adjust(left=0.10,right=0.68,top=0.88,bottom=0.12)
    x = np.arange(nks)
    for statetype in statetypes:
        col   = colors[statetype]
        ave_y = avemat[statetype].values.astype(float)
        mod_y = modmat[statetype].values.astype(float)
        deg_y = degmat[statetype].values.astype(float)
        ax.plot(x, ave_y, color=col, linestyle='-',  linewidth=2, marker='o', markersize=4)
        ax.plot(x, mod_y, color=col, linestyle='--', linewidth=2, marker='o', markersize=4)
        ax.plot(x, deg_y, color=col, linestyle=':', linewidth=2, marker='o', markersize=4)
    ax.set_xticks(x)
    ax.set_xticklabels(ks, fontsize=8)
    ax.set_xlabel('k')
    ax.set_ylabel('R²')
    ax.set_xlim(-0.3, nks - 0.7)
    color_handles = [
        mlines.Line2D([], [], color=colors[st], linewidth=2,
                      label=statetype_labels[st])
        for st in statetypes
    ]
    style_handles = [
        mlines.Line2D([], [], color='grey', linestyle='-',  linewidth=2,
                      label=ctrltype_labels['ave']),
        mlines.Line2D([], [], color='grey', linestyle='--', linewidth=2,
                      label=ctrltype_labels['mod']),
        mlines.Line2D([], [], color='grey', linestyle=':', linewidth=2,
                      label=ctrltype_labels['abs_deg'])
    ]
    fig.legend(handles=color_handles, title='Model',
               loc='upper left',
               bbox_to_anchor=(0.70, 0.88),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    fig.legend(handles=style_handles, title='Metric',
               loc='upper left',
               bbox_to_anchor=(0.70, 0.48),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    # fig.suptitle('K Comparison', fontsize=12)
    plt.savefig((outpath+'/k_SC_vs_dFC.png'),dpi=720,bbox_inches='tight')
    plt.close()

    # Cognitive version compare.
    inkey = ('/rep_testacc')
    sctype = 'fpt'
    k = '6'
    ctrltypes = ['ave','mod','abs_deg'] 
    thresval = '50'
    threstype = 'groupconsist'
    statetypes = ['SC','dFCcat','SC_dFCcat']
    nstatetypes = len(statetypes)
    septypes = ['comCFAng','comEFAng','comPCAng']
    cog_septypes = ['gCFA','gEFA','gPCA']
    ncog = len(cog_septypes)
    nrows = ncog
    ncols = nstatetypes

    # Read.
    avemat = pd.DataFrame(np.zeros((nrows,ncols)),index=cog_septypes,columns=statetypes)
    modmat = pd.DataFrame(np.zeros((nrows,ncols)),index=cog_septypes,columns=statetypes)
    degmat = pd.DataFrame(np.zeros((nrows,ncols)),index=cog_septypes,columns=statetypes)
    for cogidx in range(ncog):
        septype = septypes[cogidx]
        ccog = cog_septypes[cogidx]
        for stidx in range(nstatetypes):
            statetype = statetypes[stidx]
            for ctrltype, outmat in [('ave', avemat), ('mod', modmat), ('abs_deg', degmat)]:
                if statetype=='SC':
                    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                            + '/' + subgroup + '/'
                            '6/SC_dFC/' + sc_subgroup + '/collect/'
                            + threstype + '/' + thresval + '/KRRXF/'
                            + ctrltype + '_' + statetype + '_' + septype + '_'
                            + nrep + '_' + inner_k + '_'
                            + outer_k + '_' + sctype)
                else:
                    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                        + subgroup + '/'
                        + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                        + threstype + '/' + thresval + '/KRRXF/'
                        + ctrltype + '_' + statetype + '_' + septype + '_'
                        + '_' + nrep + '_' + inner_k + '_'
                        + outer_k + '_' + sctype)
                infile = (inpath + '/score_collect.h5')
                store = pd.HDFStore(infile, 'r')
                testacc = store.select(inkey)
                store.close()
                outmat.loc[ccog,statetype] = testacc.mean(axis=0).loc[ccog]
    
    # Plot.
    cog_labels_out = ['g (CFA)','g (EFA)','g (PCA)']
    statetype_labels = {
        'SC':        'SC',
        'dFCcat':    'dFC',
        'SC_dFCcat': 'SC & dFC',
    }
    ctrltype_labels = {
        'ave': 'AC',
        'mod': 'MC',
        'abs_deg': 'S'
    }
    colors = {
        'SC':        '#1f77b4',
        'dFCcat':    '#d62728',
        'SC_dFCcat': '#2ca02c',
    }
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.subplots_adjust(left=0.10, right=0.68, top=0.88, bottom=0.12)
    x = np.arange(ncog)                          
    for statetype in statetypes:
        col   = colors[statetype]
        ave_y = avemat[statetype].values.astype(float)
        mod_y = modmat[statetype].values.astype(float)
        deg_y = degmat[statetype].values.astype(float)
        ax.plot(x, ave_y, color=col, linestyle='-',  linewidth=2, marker='o', markersize=4)
        ax.plot(x, mod_y, color=col, linestyle='--', linewidth=2, marker='o', markersize=4)
        ax.plot(x, deg_y, color=col, linestyle=':', linewidth=2, marker='o', markersize=4)
    ax.set_xticks(x)
    ax.set_xticklabels(cog_labels_out, fontsize=8)  
    ax.set_xlabel('Cognitive Version')
    ax.set_ylabel('R²')
    ax.set_xlim(-0.3, ncog - 0.7)                
    color_handles = [
        mlines.Line2D([], [], color=colors[st], linewidth=2,
                      label=statetype_labels[st])
        for st in statetypes
    ]
    style_handles = [
        mlines.Line2D([], [], color='grey', linestyle='-',  linewidth=2,
                      label=ctrltype_labels['ave']),
        mlines.Line2D([], [], color='grey', linestyle='--', linewidth=2,
                      label=ctrltype_labels['mod']),
        mlines.Line2D([], [], color='grey', linestyle=':', linewidth=2,
                      label=ctrltype_labels['abs_deg'])              
    ]
    fig.legend(handles=color_handles, title='Model',
               loc='upper left',
               bbox_to_anchor=(0.70, 0.88),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    fig.legend(handles=style_handles, title='Metric',
               loc='upper left',
               bbox_to_anchor=(0.70, 0.48),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    plt.savefig((outpath+'/altcog_SC_vs_dFC.png'), dpi=720, bbox_inches='tight')
    plt.close()

    # SC normalization compare.
    inkey = ('/rep_testacc')
    sctypes = ['fpt','raw']
    nsctypes = len(sctypes)
    k = '6'
    ctrltypes = ['ave','mod','abs_deg'] 
    thresval = '50'
    threstype = 'groupconsist'
    statetypes = ['SC','dFCcat','SC_dFCcat']
    nstatetypes = len(statetypes)
    septype = 'comCFAng'
    cog_septype = 'gCFA'
    nrows = nsctypes
    ncols = nstatetypes

    # Read.
    avemat = pd.DataFrame(np.zeros((nrows,ncols)),index=sctypes,columns=statetypes)
    modmat = pd.DataFrame(np.zeros((nrows,ncols)),index=sctypes,columns=statetypes)
    degmat = pd.DataFrame(np.zeros((nrows,ncols)),index=sctypes,columns=statetypes)
    for scidx in range(nsctypes):
        sctype = sctypes[scidx]
        for stidx in range(nstatetypes):
            statetype = statetypes[stidx]
            for ctrltype, outmat in [('ave', avemat), ('mod', modmat), ('abs_deg', degmat)]:
                if statetype=='SC':
                    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                            + subgroup + '/'
                            '6/SC_dFC/' + sc_subgroup + '/collect/'
                            + threstype + '/' + thresval + '/KRRXFS/'
                            + ctrltype + '_' + statetype + '_' + septype + '_'
                            + '_' + nrep + '_' + inner_k + '_'
                            + outer_k + '_' + sctype)
                else:
                    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                        +  '/' + subgroup + '/'
                        + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                        + threstype + '/' + thresval + '/KRRXFS/'
                        + ctrltype + '_' + statetype + '_' + septype + '_'
                        + '_' + nrep + '_' + inner_k + '_'
                        + outer_k + '_' + sctype)
                if statetype=='dFCcat':
                    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                        + subgroup + '/'
                        + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                        + threstype + '/' + thresval + '/KRRXFS/'
                        + ctrltype + '_' + statetype + '_' + septype + '_'
                        + '_' + nrep + '_' + inner_k + '_'
                        + outer_k + '_fpt')
                infile = (inpath + '/score_collect.h5')
                store = pd.HDFStore(infile, 'r')
                testacc = store.select(inkey)
                store.close()
                outmat.loc[sctype,statetype] = testacc.mean(axis=0).loc[cog_septype]
    
    # Plot.
    sctypes_out = ['FPT','RAW']
    statetype_labels = {
        'SC': 'SC',
        'dFCcat': 'dFC',
        'SC_dFCcat': 'SC & dFC',
    }
    ctrltype_labels = {
        'ave': 'AC',
        'mod': 'MC',
        'abs_deg': 'S'
    }
    colors = {
        'SC': '#1f77b4',
        'dFCcat': '#d62728',
        'SC_dFCcat': '#2ca02c',
    }
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.subplots_adjust(left=0.10, right=0.68, top=0.88, bottom=0.12)
    x = np.arange(nsctypes)                    
    for statetype in statetypes:
        col   = colors[statetype]
        ave_y = avemat[statetype].values.astype(float)
        mod_y = modmat[statetype].values.astype(float)
        deg_y = degmat[statetype].values.astype(float)
        ax.plot(x, ave_y, color=col, linestyle='-',  linewidth=2, marker='o', markersize=4)
        ax.plot(x, mod_y, color=col, linestyle='--', linewidth=2, marker='o', markersize=4)
        ax.plot(x, deg_y, color=col, linestyle=':', linewidth=2, marker='o', markersize=4)
    ax.set_xticks(x)
    ax.set_xticklabels(sctypes_out, fontsize=8)    
    ax.set_xlabel('SC Normalization')
    ax.set_ylabel('R²')
    ax.set_xlim(-0.3, nsctypes - 0.7)          
    color_handles = [
        mlines.Line2D([], [], color=colors[st], linewidth=2,
                      label=statetype_labels[st])
        for st in statetypes
    ]
    style_handles = [
        mlines.Line2D([], [], color='grey', linestyle='-',  linewidth=2,
                      label=ctrltype_labels['ave']),
        mlines.Line2D([], [], color='grey', linestyle='--', linewidth=2,
                      label=ctrltype_labels['mod']),
        mlines.Line2D([], [], color='grey', linestyle=':', linewidth=2,
                      label=ctrltype_labels['abs_deg'])              
    ]
    fig.legend(handles=color_handles, title='Model',
               loc='upper left',
               bbox_to_anchor=(0.70, 0.88),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    fig.legend(handles=style_handles, title='Metric',
               loc='upper left',
               bbox_to_anchor=(0.70, 0.48),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    plt.savefig((outpath+'/scnorm_SC_vs_dFC.png'), dpi=720, bbox_inches='tight')
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
    statetype = 'SC_dFCcat'
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
                    +subgroup + '/'
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
    # One color per cognitive domain.
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

    # k compare.
    inkey = ('/rep_testacc')
    sctype = 'fpt'
    ks = [str(x) for x in range(2,13)] 
    nks = len(ks)
    ctrltypes = ['ave','mod','abs_deg'] 
    nctrltypes = len(ctrltypes)
    thresval = '50'
    threstype = 'groupconsist'
    statetype = 'SC_dFCcat'
    septype = 'comCFAng'
    cog_septype = ['gCFA','P24_CR','PV'] 
    ncog = len(cog_septype)
    nrows = nks
    ncols = ncog

    # Read.
    avemat = pd.DataFrame(np.zeros((nrows,ncols)),index=ks,columns=cog_septype)
    modmat = pd.DataFrame(np.zeros((nrows,ncols)),index=ks,columns=cog_septype)
    degmat = pd.DataFrame(np.zeros((nrows,ncols)),index=ks,columns=cog_septype)
    for kidx in range(nks):
        k = ks[kidx]
        for cogidx in range(ncog):
            ccog = cog_septype[cogidx]
            for ctrltype, outmat in [('ave', avemat), ('mod', modmat), ('abs_deg', degmat)]:
                inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                    + subgroup + '/'
                    + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                    + threstype + '/' + thresval + '/KRRXF/'
                    + ctrltype + '_' + statetype + '_' + septype + '_'
                    + nrep + '_' + inner_k + '_'
                    + outer_k + '_' + sctype)
                infile = (inpath + '/score_collect.h5')
                store = pd.HDFStore(infile, 'r')
                testacc = store.select(inkey)
                store.close()
                outmat.loc[k,ccog] = testacc.mean(axis=0).loc[ccog]

    # Plot.
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
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.subplots_adjust(left=0.10, right=0.68, top=0.88, bottom=0.12)
    x = np.arange(nks)
    for ccog in cog_septype:
        col   = cog_colors[ccog]
        ave_y = avemat[ccog].values.astype(float)
        mod_y = modmat[ccog].values.astype(float)
        deg_y = degmat[ccog].values.astype(float)
        ax.plot(x, ave_y, color=col, linestyle='-',  linewidth=2, marker='o', markersize=4)
        ax.plot(x, mod_y, color=col, linestyle='--', linewidth=2, marker='o', markersize=4)
        ax.plot(x, deg_y, color=col, linestyle=':', linewidth=2, marker='o', markersize=4)
    ax.set_xticks(x)
    ax.set_xticklabels(ks, fontsize=8)
    ax.set_xlabel('k')
    ax.set_ylabel('R²')
    ax.set_xlim(-0.3, nks - 0.7)
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
               bbox_to_anchor=(0.70, 0.88),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    fig.legend(handles=style_handles, title='Metric',
               loc='upper left',
               bbox_to_anchor=(0.70, 0.38),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    plt.savefig((outpath+'/k_g_vs_gF_vs_gC.png'),dpi=720,bbox_inches='tight')
    plt.close()

    # Cognitive version compare.
    inkey = ('/rep_testacc')
    sctype = 'fpt'
    k = '6'
    ctrltypes = ['ave','mod','abs_deg'] 
    nctrltypes = len(ctrltypes)
    thresval = '50'
    threstype = 'groupconsist'
    statetype = 'SC_dFCcat'
    septypes = ['comCFAng','comEFAng','comPCAng','comNT']
    nseptypes = len(septypes)
    cog_septypes = [['gCFA','P24_CR','PV'],
                    ['gEFA'],
                    ['gPCA2'],
                    ['NF','NC']]
    cog_comtypes = ['gCFA','P24_CR','PV','gEFA','gPCA','NF','NC'] 
    ncog = len(cog_comtypes)
    cog_labels_out = ['g (CFA)','gF (PMAT24)','gC (PicVoc)',
                      'g (EFA)','g (PCA)','gF (NIH)','gC (NIH)']
    nrows = nctrltypes
    ncols = ncog

    # Read.
    allmat = pd.DataFrame(np.zeros((nrows,ncols)),index=ctrltypes,columns=cog_comtypes)
    for seidx in range(nseptypes):
        septype = septypes[seidx]
        cseptypes = cog_septypes[seidx]
        for ccog in cseptypes:
            for ctrltype in ctrltypes:
                inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                    + '/' + subgroup + '/'
                    + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                    + threstype + '/' + thresval + '/KRRXFS/'
                    + ctrltype + '_' + statetype + '_' + septype + '_'
                    + '_' + nrep + '_' + inner_k + '_'
                    + outer_k + '_' + sctype)
                infile = (inpath + '/score_collect.h5')
                store = pd.HDFStore(infile, 'r')
                testacc = store.select(inkey)
                store.close()
                allmat.loc[ctrltype,ccog] = testacc.mean(axis=0).loc[ccog]
    
    # Plot.
    ctrltype_labels = {
        'ave': 'AC',
        'mod': 'MC',
        'abs_deg': 'S'
    }
    colors = {
        'ave': '#1f77b4',
        'mod': '#d62728',
        'abs_deg':  '#2ca02c',
    }
    group_boundaries = np.cumsum([len(g) for g in cog_septypes[:-1]]) - 0.5
    group_names = ['CFA','EFA','PCA','NIH']
    group_sizes = [len(g) for g in cog_septypes]
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.subplots_adjust(left=0.08, right=0.68, top=0.88, bottom=0.18)
    x = np.arange(ncog)
    for ctrltype in ctrltypes:
        col = colors[ctrltype]
        y   = allmat.loc[ctrltype].values.astype(float)
        if ctrltype=='ave':
            cline = '-'
        elif ctrltype=='mod':
            cline = '--'
        elif ctrltype=='abs_deg':
            cline = ':'
        ax.plot(x, y, color=col,
                linestyle=cline,
                linewidth=2, marker='o', markersize=4,
                label=ctrltype_labels[ctrltype])
    for boundary in group_boundaries:
        ax.axvline(boundary, color='grey', linestyle=':', linewidth=0.8, zorder=0)
    start = 0
    for gname, gsize in zip(group_names, group_sizes):
        mid = start + (gsize - 1) / 2
        ax.text(mid, -0.13, gname, ha='center', va='top',
                transform=ax.get_xaxis_transform(), fontsize=8, color='grey')
        start += gsize
    ax.set_xticks(x)
    ax.set_xticklabels(cog_labels_out, fontsize=7, rotation=30, ha='right')
    ax.set_xlabel('Cognitive Version',labelpad=10)
    ax.set_ylabel('R²')
    ax.set_xlim(-0.3, ncog - 0.7)
    style_handles = []
    for ct in ctrltypes:
        if ctrltype=='ave':
            cline = '-'
        elif ctrltype=='mod':
            cline = '--'
        elif ctrltype=='abs_deg':
            cline = ':'
        style_handles.append(mlines.Line2D([], [], color=colors[ct], linestyle=cline,
                            linewidth=2, label=ctrltype_labels[ct]))
    fig.legend(handles=style_handles, title='Metric',
               loc='upper left',
               bbox_to_anchor=(0.70, 0.88),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    plt.savefig((outpath+'/altcog_g_vs_gF_vs_gC.png'), dpi=720, bbox_inches='tight')
    plt.close()

    # SC normalization compare.
    inkey = ('/rep_testacc')
    sctypes = ['fpt','raw']
    nsctypes = len(sctypes)
    k = '6'
    ctrltypes = ['ave','mod','abs_deg'] 
    thresval = '50'
    threstype = 'groupconsist'
    statetype = 'SC_dFCcat'
    septype = 'comCFAng'
    cog_septype = ['gCFA','P24_CR','PV']
    ncog = len(cog_septype)
    nrows = nsctypes
    ncols = ncog

    # Read.
    avemat = pd.DataFrame(np.zeros((nrows, ncols)), index=sctypes, columns=cog_septype)
    modmat = pd.DataFrame(np.zeros((nrows, ncols)), index=sctypes, columns=cog_septype)
    degmat = pd.DataFrame(np.zeros((nrows, ncols)), index=sctypes, columns=cog_septype)
    for scidx in range(nsctypes):
        sctype = sctypes[scidx]
        for cogidx in range(ncog):
            ccog = cog_septype[cogidx]
            for ctrltype, outmat in [('ave', avemat), ('mod', modmat), ('abs_deg', degmat)]:  
                inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                    + '/' + subgroup + '/'
                    + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                    + threstype + '/' + thresval + '/KRRXFS/'
                    + ctrltype + '_' + statetype + '_' + septype + '_'
                    + nrep + '_' + inner_k + '_'
                    + outer_k + '_' + sctype)
                infile = (inpath + '/score_collect.h5')
                store = pd.HDFStore(infile, 'r')
                testacc = store.select(inkey)
                store.close()
                outmat.loc[sctype, ccog] = testacc.mean(axis=0).loc[ccog]
    
    # Plot.
    sctypes_out = ['FPT','RAW']
    cog_labels = {
        'gCFA':    'g',
        'P24_CR':  'gF',
        'PV':      'gC'
    }
    ctrltype_labels = {
        'ave': 'AC',
        'mod': 'MC',
        'abs_deg': 'S'
    }
    cog_colors = {
        'gCFA':    '#1f77b4',
        'P24_CR':  '#d62728',
        'PV':      '#2ca02c'
    }
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.subplots_adjust(left=0.10, right=0.68, top=0.88, bottom=0.12)
    x = np.arange(nsctypes)
    for ccog in cog_septype:
        col   = cog_colors[ccog]
        ave_y = avemat[ccog].values.astype(float)
        mod_y = modmat[ccog].values.astype(float)
        deg_y = degmat[ccog].values.astype(float)
        ax.plot(x, ave_y, color=col, linestyle='-',  linewidth=2, marker='o', markersize=4)
        ax.plot(x, mod_y, color=col, linestyle='--', linewidth=2, marker='o', markersize=4)
        ax.plot(x, deg_y, color=col, linestyle=':', linewidth=2, marker='o', markersize=4)
    ax.set_xticks(x)
    ax.set_xticklabels(sctypes_out, fontsize=8)
    ax.set_xlabel('SC Normalization')
    ax.set_ylabel('R²')
    ax.set_xlim(-0.3, nsctypes - 0.7)
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
               bbox_to_anchor=(0.70, 0.88),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    fig.legend(handles=style_handles, title='Metric',
               loc='upper left',
               bbox_to_anchor=(0.70, 0.38),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    plt.savefig((outpath+'/scnorm_g_vs_gF_vs_gC.png'), dpi=720, bbox_inches='tight')
    plt.close()

    # ---------------------------- Principal Gradient ---------------------------- #

    # Read in principal gradient.
    nroi = 360
    infile = ('../outputs/r_sFC/dr_full/none/0/sFC_gradients.csv')
    sFCgr_all = pd.read_csv(infile,header=None)
    infile = ('../outputs/r_sFC/dr_full/none/0/sFC_gradients_flip.csv')
    sFCgr_flip = pd.read_csv(infile,header=None).values.tolist()[0]
    ngr = len(sFCgr_flip)
    for gidx in range(ngr):
        if sFCgr_flip[gidx] == 'T':
            sFCgr_all.iloc[:,gidx] = -sFCgr_all.iloc[:,gidx]
    sFCgr_all.index = [('r'+str(ridx+1)) for ridx in range(nroi)]
    sFCgr = sFCgr_all.iloc[:,0]

    # Threshold compare.
    inkey = ('/rep_covha_feat_gCFA')
    sctype = 'fpt'
    k = '6' 
    nk = int(k)
    klabs = [f's{str(x+1)}' for x in range(nk)] # One line for each state.
    klabs_out = [f'S{str(x+1)}' for x in range(nk)] # What to display.
    ctrltypes = ['ave','mod','abs_deg']
    nctrltypes = len(ctrltypes)
    thresvals = ['10','20','30','40','50','60','70','80','90'] 
    nthresvals = len(thresvals)
    consist_labs = ['none_0'] + [('groupconsist_'+x) for x in thresvals]
    str_labs = ['none_0'] + [('groupstr_'+x) for x in thresvals]
    threslabs = consist_labs + [('groupstr_'+x) for x in thresvals]
    nthres = len(threslabs)
    statetype = 'SC_dFCcat'
    septype = 'comCFAng'
    nrows = nthres
    ncols = nk

    # Read, get each dFC state and find correlation with gradient.
    avemat = pd.DataFrame(np.zeros((nrows,ncols)),index=threslabs,columns=klabs)
    modmat = pd.DataFrame(np.zeros((nrows,ncols)),index=threslabs,columns=klabs)
    degmat = pd.DataFrame(np.zeros((nrows,ncols)),index=threslabs,columns=klabs)
    for thidx in range(nthres):
        cthres = threslabs[thidx]
        threstype,thresval = cthres.split('_')
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
            featimp = store.select(inkey)
            store.close()
            allfeat = featimp.mean(axis=0)
            for kidx in range(nk):
                cstate = klabs[kidx]
                fcols = [x for x in allfeat.index if f'{cstate}_' in x]
                cstate_feat = allfeat[fcols]
                gcols = [x.replace(f'{cstate}_','') for x in fcols]
                cgr = sFCgr[gcols]
                outmat.loc[cthres,cstate] = spearmanr(cstate_feat,cgr).statistic
    
    # Plot.
    ax_titles = {
        'groupconsist': 'Consistency',
        'groupstr':     'Strength',
    }
    ctrltype_labels = {
        'ave': 'AC',
        'mod': 'MC',
        'abs_deg': 'S'
    }
    cmap = plt.cm.get_cmap('tab10',nk)
    state_colors = {klabs[i]: cmap(i) for i in range(nk)}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    fig.subplots_adjust(wspace=0.08, right=0.78)
    x_shared = np.arange(len(thresvals) + 1)
    for ax, group_labs, group_key in [(ax1, consist_labs, 'groupconsist'),
                                      (ax2, str_labs,     'groupstr')]:
        for kidx in range(nk):
            cstate = klabs[kidx]
            col   = state_colors[cstate]
            ave_y = np.array([avemat.loc[tl, cstate] for tl in group_labs], dtype=float)
            mod_y = np.array([modmat.loc[tl, cstate] for tl in group_labs], dtype=float)
            deg_y = np.array([degmat.loc[tl, cstate] for tl in group_labs], dtype=float)
            ax.plot(x_shared, ave_y, color=col, linestyle='-',  linewidth=2, marker='o', markersize=4)
            ax.plot(x_shared, mod_y, color=col, linestyle='--', linewidth=2, marker='o', markersize=4)
            ax.plot(x_shared, deg_y, color=col, linestyle=':', linewidth=2, marker='o', markersize=4)
        ax.axhline(0, color='black', linewidth=0.8, linestyle=':', zorder=0)
        tick_labels = ['0'] + thresvals
        ax.set_xticks(x_shared)
        ax.set_xticklabels(tick_labels, fontsize=8)
        ax.set_xlabel('Threshold (% Removed)')
        ax.set_title(ax_titles[group_key], fontsize=10)
        ax.set_xlim(-0.3, len(thresvals) + 0.3)
    ax1.set_ylabel('Principal Gradient RHO')
    state_handles = [
        mlines.Line2D([], [], color=state_colors[klabs[i]], linewidth=2,
                      label=klabs_out[i])
        for i in range(nk)
    ]
    style_handles = [
        mlines.Line2D([], [], color='grey', linestyle='-',  linewidth=2,
                      label=ctrltype_labels['ave']),
        mlines.Line2D([], [], color='grey', linestyle='--', linewidth=2,
                      label=ctrltype_labels['mod']),
        mlines.Line2D([], [], color='grey', linestyle=':', linewidth=2,
                      label=ctrltype_labels['abs_deg'])              
    ]
    fig.legend(handles=state_handles, title='State',
               loc='upper left',
               bbox_to_anchor=(0.80, 0.92),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    fig.legend(handles=style_handles, title='Metric',
               loc='upper left',
               bbox_to_anchor=(0.80, 0.30),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    plt.savefig((outpath+'/threshold_sFCgr.png'), dpi=720, bbox_inches='tight')
    plt.close()

    # k compare.
    inkey = ('/rep_covha_feat_gCFA')
    sctype = 'fpt'
    ks = [str(x) for x in range(2,13)] 
    nks = len(ks)
    maxk = 12
    klabs = [f's{str(x+1)}' for x in range(maxk)] 
    klabs_out = [f'S{str(x+1)}' for x in range(maxk)] 
    ctrltypes = ['ave','mod','abs_deg'] 
    nctrltypes = len(ctrltypes)
    thresval = '50'
    threstype = 'groupconsist'
    statetype = 'SC_dFCcat'
    septype = 'comCFAng'
    nrows = nks
    ncols = maxk

    # Read, get each dFC state and find correlation with gradient.
    avemat = pd.DataFrame(np.full((nrows,ncols),np.nan),index=ks,columns=klabs)
    modmat = pd.DataFrame(np.full((nrows,ncols),np.nan),index=ks,columns=klabs)
    degmat = pd.DataFrame(np.full((nrows,ncols),np.nan),index=ks,columns=klabs)
    for kidx in range(nks):
        k = ks[kidx]
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
            featimp = store.select(inkey)
            store.close()
            if int(k) >= 11:
                allfeat = featimp.T.mean(axis=0)
            else:
                allfeat = featimp.mean(axis=0)
            for stateidx in range(int(k)):
                cstate = klabs[stateidx]
                fcols = [x for x in allfeat.index if f'{cstate}_' in x]
                cstate_feat = allfeat[fcols]
                gcols = [x.replace(f'{cstate}_','') for x in fcols]
                cgr = sFCgr[gcols]
                outmat.loc[k,cstate] = spearmanr(cstate_feat,cgr).statistic

    # Plot.
    ctrltype_labels = {
        'ave': 'AC',
        'mod': 'MC',
        'abs_deg': 'S',
    }
    def saturate(color, saturation=0.9):
        r, g, b, *_ = color
        h, l, s = colorsys.rgb_to_hls(r, g, b)
        r2, g2, b2 = colorsys.hls_to_rgb(h, l, saturation)
        return (r2, g2, b2)
    raw_cmap = plt.cm.get_cmap('Set3', maxk)
    state_colors = {klabs[i]: saturate(raw_cmap(i)) for i in range(maxk)}
    state_colors['s2'] = '#e6b800'  # Replace faint yellow with a deeper gold.
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.subplots_adjust(left=0.10, right=0.68, top=0.88, bottom=0.12)
    x = np.arange(nks)
    for i in range(maxk):
        cstate = klabs[i]
        col   = state_colors[cstate]
        ave_y = avemat[cstate].values.astype(float)  # NaN where state doesn't appear yet.
        mod_y = modmat[cstate].values.astype(float)
        deg_y = degmat[cstate].values.astype(float)
        ax.plot(x, ave_y, color=col, linestyle='-',  linewidth=2, marker='o', markersize=4)
        ax.plot(x, mod_y, color=col, linestyle='--', linewidth=2, marker='o', markersize=4)
        ax.plot(x, deg_y, color=col, linestyle=':', linewidth=2, marker='o', markersize=4)
    ax.axhline(0, color='black', linewidth=0.8, linestyle=':', zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels(ks, fontsize=8)
    ax.set_xlabel('k')
    ax.set_ylabel('Principal Gradient RHO')
    ax.set_xlim(-0.3, nks - 0.7)
    state_handles = [
        mlines.Line2D([], [], color=state_colors[klabs[i]], linewidth=2,
                      label=klabs_out[i])
        for i in range(maxk)
    ]
    style_handles = [
        mlines.Line2D([], [], color='grey', linestyle='-',  linewidth=2,
                      label=ctrltype_labels['ave']),
        mlines.Line2D([], [], color='grey', linestyle='--', linewidth=2,
                      label=ctrltype_labels['mod']),
        mlines.Line2D([], [], color='grey', linestyle=':', linewidth=2,
                      label=ctrltype_labels['abs_deg'])
    ]
    fig.legend(handles=state_handles, title='State',
               loc='upper left',
               bbox_to_anchor=(0.70, 0.88),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    fig.legend(handles=style_handles, title='Metric',
               loc='upper left',
               bbox_to_anchor=(0.70, 0.18),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    plt.savefig((outpath+'/k_sFCgr.png'), dpi=720, bbox_inches='tight')
    plt.close()

    # Cognitive version compare.
    sctype = 'fpt'
    k = '6' 
    nk = int(k)
    klabs = [f's{str(x+1)}' for x in range(nk)] 
    klabs_out = [f'S{str(x+1)}' for x in range(nk)] 
    ctrltypes = ['ave','mod','abs_deg'] 
    nctrltypes = len(ctrltypes)
    thresval = '50'
    threstype = 'groupconsist'
    statetype = 'SC_dFCcat'
    septypes = ['comCFAng','comEFAng','comPCAng']
    cog_septypes = ['gCFA','gEFA','gPCA']
    ncog = len(cog_septypes)
    nrows = ncog
    ncols = nk

    # Read, get each dFC state and find correlation with gradient.
    avemat = pd.DataFrame(np.zeros((nrows,ncols)),index=cog_septypes,columns=klabs)
    modmat = pd.DataFrame(np.zeros((nrows,ncols)),index=cog_septypes,columns=klabs)
    degmat = pd.DataFrame(np.zeros((nrows,ncols)),index=cog_septypes,columns=klabs)
    for cogidx in range(ncog):
        septype = septypes[cogidx]
        ccog = cog_septypes[cogidx]
        for ctrltype, outmat in [('ave', avemat), ('mod', modmat), ('abs_deg', degmat)]:
            inkey = ('/rep_covha_feat_'+ccog)
            inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                + subgroup + '/'
                + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                + threstype + '/' + thresval + '/KRRXFS/'
                + ctrltype + '_' + statetype + '_' + septype + '_'
                + nrep + '_' + inner_k + '_'
                + outer_k + '_' + sctype)
            infile = (inpath + '/score_collect.h5')
            store = pd.HDFStore(infile, 'r')
            featimp = store.select(inkey)
            store.close()
            allfeat = featimp.mean(axis=0)
            for kidx in range(nk):
                cstate = klabs[kidx]
                fcols = [x for x in allfeat.index if f'{cstate}_' in x]
                cstate_feat = allfeat[fcols]
                gcols = [x.replace(f'{cstate}_','') for x in fcols]
                cgr = sFCgr[gcols]
                outmat.loc[ccog,cstate] = spearmanr(cstate_feat,cgr).statistic
    
    # Plot.
    cog_labels_out = ['g (CFA)','g (EFA)','g (PCA)'] 
    ctrltype_labels = {
        'ave': 'AC',
        'mod': 'MC',
        'abs_deg': 'S'
    }
    state_colors = {klabs[i]: plt.cm.tab10(i) for i in range(nk)}
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.subplots_adjust(left=0.10, right=0.68, top=0.88, bottom=0.12)
    x = np.arange(ncog)                              
    for kidx in range(nk):
        cstate = klabs[kidx]
        col   = state_colors[cstate]
        ave_y = avemat.loc[cog_septypes, cstate].values.astype(float)
        mod_y = modmat.loc[cog_septypes, cstate].values.astype(float)
        deg_y = degmat.loc[cog_septypes, cstate].values.astype(float)
        ax.plot(x, ave_y, color=col, linestyle='-',  linewidth=2, marker='o', markersize=4)
        ax.plot(x, mod_y, color=col, linestyle='--', linewidth=2, marker='o', markersize=4)
        ax.plot(x, deg_y, color=col, linestyle=':', linewidth=2, marker='o', markersize=4)
    ax.axhline(0, color='black', linewidth=0.8, linestyle=':', zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels(cog_labels_out, fontsize=8)
    ax.set_xlabel('Cognitive Version')
    ax.set_ylabel('Principal Gradient RHO')
    ax.set_xlim(-0.3, ncog - 0.7)
    state_handles = [
        mlines.Line2D([], [], color=state_colors[klabs[i]], linewidth=2,
                      label=klabs_out[i])
        for i in range(nk)
    ]
    style_handles = [
        mlines.Line2D([], [], color='grey', linestyle='-',  linewidth=2,
                      label=ctrltype_labels['ave']),
        mlines.Line2D([], [], color='grey', linestyle='--', linewidth=2,
                      label=ctrltype_labels['mod']),
        mlines.Line2D([], [], color='grey', linestyle=':', linewidth=2,
                      label=ctrltype_labels['abs_deg'])
    ]
    fig.legend(handles=state_handles, title='State',
               loc='upper left',
               bbox_to_anchor=(0.70, 0.88),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    fig.legend(handles=style_handles, title='Metric',
               loc='upper left',
               bbox_to_anchor=(0.70, 0.30),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    plt.savefig((outpath+'/altcog_sFCgr.png'), dpi=720, bbox_inches='tight')
    plt.close()

    # SC normalizaton compare.
    sctypes = ['fpt','raw']
    nsctypes = len(sctypes)
    k = '6' 
    nk = int(k)
    klabs = [f's{str(x+1)}' for x in range(nk)] 
    klabs_out = [f'S{str(x+1)}' for x in range(nk)] 
    ctrltypes = ['ave','mod','abs_deg']
    nctrltypes = len(ctrltypes)
    thresval = '50'
    threstype = 'groupconsist'
    statetype = 'SC_dFCcat'
    septype = 'comCFAng'
    cog_septype = 'gCFA'
    nrows = nsctypes
    ncols = nk

    # Read, get each dFC state and find correlation with gradient.
    avemat = pd.DataFrame(np.zeros((nrows,ncols)),index=sctypes,columns=klabs)
    modmat = pd.DataFrame(np.zeros((nrows,ncols)),index=sctypes,columns=klabs)
    degmat = pd.DataFrame(np.zeros((nrows,ncols)),index=sctypes,columns=klabs)
    for scidx in range(nsctypes):
        sctype = sctypes[scidx]
        for ctrltype, outmat in [('ave', avemat), ('mod', modmat), ('abs_deg', degmat)]:
            inkey = ('/rep_covha_feat_gCFA')
            inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                + subgroup + '/'
                + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                + threstype + '/' + thresval + '/KRRXFS/'
                + ctrltype + '_' + statetype + '_' + septype + '_'
                + nrep + '_' + inner_k + '_'
                + outer_k + '_' + sctype)
            infile = (inpath + '/score_collect.h5')
            store = pd.HDFStore(infile, 'r')
            featimp = store.select(inkey)
            store.close()
            allfeat = featimp.mean(axis=0)
            for kidx in range(nk):
                cstate = klabs[kidx]
                fcols = [x for x in allfeat.index if f'{cstate}_' in x]
                cstate_feat = allfeat[fcols]
                gcols = [x.replace(f'{cstate}_','') for x in fcols]
                cgr = sFCgr[gcols]
                outmat.loc[sctype,cstate] = spearmanr(cstate_feat,cgr).statistic
    
    # Plot.
    sctypes_out = ['FPT','RAW'] 
    ctrltype_labels = {
        'ave': 'AC',
        'mod': 'MC',
        'abs_deg': 'S'
    }
    state_colors = {klabs[i]: plt.cm.tab10(i) for i in range(nk)}
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.subplots_adjust(left=0.10, right=0.68, top=0.88, bottom=0.12)
    x = np.arange(nsctypes)                          
    for kidx in range(nk):
        cstate = klabs[kidx]
        col   = state_colors[cstate]
        ave_y = avemat.loc[sctypes, cstate].values.astype(float)
        mod_y = modmat.loc[sctypes, cstate].values.astype(float)
        deg_y = degmat.loc[sctypes, cstate].values.astype(float)
        ax.plot(x, ave_y, color=col, linestyle='-',  linewidth=2, marker='o', markersize=4)
        ax.plot(x, mod_y, color=col, linestyle='--', linewidth=2, marker='o', markersize=4)
        ax.plot(x, deg_y, color=col, linestyle=':', linewidth=2, marker='o', markersize=4)
    ax.axhline(0, color='black', linewidth=0.8, linestyle=':', zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels(sctypes_out, fontsize=8)      
    ax.set_xlabel('SC Normalization')
    ax.set_ylabel('Principal Gradient RHO')
    ax.set_xlim(-0.3, nsctypes - 0.7)               
    state_handles = [
        mlines.Line2D([], [], color=state_colors[klabs[i]], linewidth=2,
                      label=klabs_out[i])
        for i in range(nk)
    ]
    style_handles = [
        mlines.Line2D([], [], color='grey', linestyle='-',  linewidth=2,
                      label=ctrltype_labels['ave']),
        mlines.Line2D([], [], color='grey', linestyle='--', linewidth=2,
                      label=ctrltype_labels['mod']),
        mlines.Line2D([], [], color='grey', linestyle=':', linewidth=2,
                      label=ctrltype_labels['abs_deg'])              
    ]
    fig.legend(handles=state_handles, title='State',
               loc='upper left',
               bbox_to_anchor=(0.70, 0.88),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    fig.legend(handles=style_handles, title='Metric',
               loc='upper left',
               bbox_to_anchor=(0.70, 0.30),
               bbox_transform=fig.transFigure,
               borderaxespad=0, framealpha=0.8)
    plt.savefig((outpath+'/scnorm_sFCgr.png'), dpi=720, bbox_inches='tight')
    plt.close()
    print('Parameter variation done.')
