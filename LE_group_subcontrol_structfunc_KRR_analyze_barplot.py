# -*- coding: utf-8 -*-
"""
For a given number of repetitions, inner loop number, and outer loop number for CV,
for the main models of interest including SC and dFC, and the supplementary models 
including sFC as well, plot bar plots with along with annotations for numerical 
characteristics including percentage contributions in the models.
Output:
vert_SC_dFC_barplots.png Bar plots for the main models including SC and dFC.
vert_SC_sFC_dFC_barplots.png Bar plots for the supplementary models including sFC as well.

Usage: 
    LE_group_subcontrol_structfunc_KRR_analyze_barplot.py <nrep> <inner_k> <outer_k> 
    
Arguments:
    
    <nrep> Number of CV repetitions
    <inner_k> Inner K in K-fold CV hyperparameter search
    <outer_k> Outer K in K-fold CV
    <nrep> Number of permutations

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
               subgroup+'/SC_dFC_comkthres/'+sc_subgroup+'/KRRXFS/'+nrep+'_'+inner_k+'_'+outer_k+'/')
    os.makedirs(outpath,exist_ok=True)

    # Set parameters corresponding to the models of interest.
    k = '6'
    sctype = 'fpt'
    threstype = 'groupconsist'
    thresval = '50'
    septype = 'comCFAng'
    coglist = ['gCFA','P24_CR','PV']
    coglist_lab = ['g','gF','gC']
    ncog = len(coglist)
    controltypes = ['ave','mod','abs_deg'] 
    controltypes_lab = ['AC','MC','S']
    nctrl = len(controltypes)

    # -------------------------------- SC and dFC -------------------------------- #

    # Read in.
    plotmn_dict_main = {}
    plotp_dict = {}
    plotperc_dict_main = {}
    statetypes_main = ['SC_dFCcat','dFCcat','SC']
    statetypes_lab_main = ['SC & dFC','dFC','SC']
    statetype_perc_main = {'SC_dFCcat':['sc','dFC']}
    statetype_perc_lab_main = {'SC_dFCcat':['SC','dFC']}
    nstatetypes = len(statetypes_main)
    for ccog in coglist:
        plotmn_main = pd.DataFrame(np.zeros((nctrl,nstatetypes)),index=controltypes,columns=statetypes_main)
        plotp = pd.DataFrame(np.zeros((nctrl,nstatetypes)),index=controltypes,columns=statetypes_main)
        for ctrltype in controltypes:
            for statetype in statetypes_main:

                # Read R2 and p-value.
                inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                                        + subgroup + '/'
                                        + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                                        + threstype + '/' + thresval + '/KRRXFS/'
                                        + ctrltype + '_' + statetype + '_' + septype + '_'
                                        + nrep + '_' + inner_k + '_'
                                        + outer_k + '_' + sctype)
                infile = (inpath+'/testacc_FDR.csv')
                inmat = pd.read_csv(infile,index_col=0)
                ckey = (f'{ctrltype} {statetype} {ccog}')
                plotmn_main.loc[ctrltype,statetype] = inmat.loc[ckey,'Accuracy']
                plotp.loc[ctrltype,statetype] = inmat.loc[ckey,'OneP_FDR']

                # Read percentages if they exist.
                if statetype in statetype_perc_main.keys():
                    infile = (inpath+'/splitacc_strict_perc.csv')
                    inmat = pd.read_csv(infile,index_col=0)
                    cperclabs = statetype_perc_main[statetype]
                    plotperc_dict_main[ckey] = inmat.loc[cperclabs,ccog]

        # Append the list.
        plotmn_dict_main[ccog] = plotmn_main
        plotp_dict[ccog] = plotp

    # ---------------------------- SC and dFC and sFC ---------------------------- #
    
    # Read in.
    plotmn_dict_supp = {}
    plotperc_dict_supp = {}
    statetypes_supp = ['SC_sFC_dFCcat','SC_sFC','sFC']
    statetypes_lab_supp = ['SC & sFC & dFC','SC & sFC','sFC']
    statetype_perc_supp = {'SC_sFC_dFCcat':['sc','sFC','dFC'],
                        'SC_sFC':['sc','sFC']}
    statetype_perc_lab_supp = {'SC_sFC_dFCcat':['SC','sFC','dFC'],
                            'SC_sFC':['SC','sFC']}
    nstatetypes = len(statetypes_supp)
    for ccog in coglist:
        plotmn = pd.DataFrame(np.zeros((nctrl,nstatetypes)),index=controltypes,columns=statetypes_supp)
        for ctrltype in controltypes:
            for statetype in statetypes_supp:

                # Read R2.
                inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                                        + subgroup + '/'
                                        + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                                        + threstype + '/' + thresval + '/KRRXF/'
                                        + ctrltype + '_' + statetype + '_' + septype + '_'
                                        + nrep + '_' + inner_k + '_'
                                        + outer_k + '_' + sctype)
                infile = (inpath+'/testacc_noP.csv')
                inmat = pd.read_csv(infile,index_col=0)
                ckey = (f'{ctrltype} {statetype} {ccog}')
                plotmn.loc[ctrltype,statetype] = inmat.loc[ccog,'Accuracy']

                # Read percentages if they exist.
                if statetype in statetype_perc_supp.keys():
                    infile = (inpath+'/splitacc_strict_perc.csv')
                    inmat = pd.read_csv(infile,index_col=0)
                    cperclabs = statetype_perc_supp[statetype]
                    plotperc_dict_supp[ckey] = inmat.loc[cperclabs,ccog]

        # Append the list.
        plotmn_dict_supp[ccog] = plotmn

    # ----------------------------------- Plot ----------------------------------- #

    # Compute a single y-axis range across both plots for comparison.
    all_vals = []
    for ccog in coglist:
        all_vals.extend(plotmn_dict_main[ccog].values.flatten())
        all_vals.extend(plotmn_dict_supp[ccog].values.flatten())
    global_ymin = min(0,np.min(all_vals))  
    global_ymax = np.max(all_vals)
    global_yrange = global_ymax - global_ymin
    shared_ylim = (global_ymin, global_ymax + 0.35 * global_yrange)

    # Plot one plot for each cognitive variable. Plot bar plots for each R2,
    # where AC, MC, and S are grouped together. SC, dFC, and SC & dFC are spaced apart.
    # Annotate on top of the bar plots from top-down - first a star for significance 
    # if it is significant, the R2 value, then the percentages if 
    # they exist (on top of each other), then the p-value. Make it < if it is the
    # p-value from FDR-correction of the the p-value where no null value was higher.
    minp = np.min([np.min(plotp_dict[x]) for x in coglist])
    nstate = nstatetypes
    ncontrol = nctrl
    gap = 0.7
    barwidth = 0.8 / ncontrol
    xpos = []
    curr = 0
    for stateidx in range(nstate):
        for ctrlidx in range(ncontrol):
            xpos.append(curr)
            curr += 1
        curr += gap
    xpos = np.array(xpos)
    ctrl_colors = {'ave': '#1f77b4', 'mod': '#d62728', 'abs_deg': '#2ca02c'}
    sig_thresh = 0.05
    fontsize_annot = 10
    fig, axes = plt.subplots(ncog, 1, figsize=(12,4.5*ncog),sharex=True,sharey=True)
    for cidx, ccog in enumerate(coglist):
        ax = axes[cidx]
        plotmn = plotmn_dict_main[ccog]
        plotp = plotp_dict[ccog]
        for ctrlidx, ctrltype in enumerate(controltypes):
            cvals = []
            barpos = []
            idx = 0
            for stateidx, statetype in enumerate(statetypes_main):
                barpos.append(xpos[idx + ctrlidx])
                cvals.append(plotmn.loc[ctrltype, statetype])
                idx += ncontrol
            ax.bar(barpos, cvals, width=barwidth,
                    label=controltypes_lab[ctrlidx],
                    color=ctrl_colors[ctrltype],
                    edgecolor='black', linewidth=0.5)
            idx = 0
            for stateidx, statetype in enumerate(statetypes_main):
                val = plotmn.loc[ctrltype, statetype]
                pval = plotp.loc[ctrltype, statetype]
                cx = xpos[idx + ctrlidx]
                idx += ncontrol
                bar_top = max(val, 0)
                x_offset = -9
                pt_offset = 3
                if pval == minp:
                    plab_here = f'p < {pval:.3f}'
                else:
                    plab_here = f'p = {pval:.3f}'
                ax.annotate(plab_here, xy=(cx,bar_top), xycoords='data',
                            xytext=(x_offset,pt_offset), textcoords='offset points',
                            ha='left', va='bottom', fontsize=fontsize_annot)
                pt_offset += 13
                ckey = f'{ctrltype} {statetype} {ccog}'
                if statetype in statetype_perc_main.keys():
                    percs = plotperc_dict_main[ckey]
                    perclabs = statetype_perc_lab_main[statetype]
                    for plabel, praw in zip(perclabs, percs.index):
                        ax.annotate(f'{plabel}: {percs[praw]:.1f}%', xy=(cx,bar_top), xycoords='data',
                                    xytext=(x_offset,pt_offset), textcoords='offset points',
                                    ha='left', va='bottom', fontsize=fontsize_annot)
                        pt_offset += 13
                ax.annotate(f'{val:.3f}', xy=(cx,bar_top), xycoords='data',
                            xytext=(x_offset,pt_offset), textcoords='offset points',
                            ha='left', va='bottom', fontsize=fontsize_annot, fontweight='bold')
                pt_offset += 15
                if pval < sig_thresh:
                    ax.annotate('*', xy=(cx,bar_top), xycoords='data',
                                xytext=(0,pt_offset), textcoords='offset points',
                                ha='center', va='bottom', fontsize=fontsize_annot*1.5,
                                fontweight='bold', color='red')
        ax.axhline(0, color='black', linewidth=0.8)
        if cidx == 1:
            ax.set_ylabel('R²')
        ax.text(-0.10, 0.5, coglist_lab[cidx],transform=ax.transAxes,
            ha='center', va='center', rotation=0, fontsize=12, fontweight='bold')
        ax.set_ylim(*shared_ylim)
        if cidx == 0:
            ax.legend(title='Metric', loc='upper right')
    block_centers = []
    idx = 0
    for stateidx in range(nstate):
        block_xpos = xpos[idx:idx+ncontrol]
        block_centers.append(block_xpos.mean())
        idx += ncontrol
    axes[-1].set_xticks(block_centers)
    axes[-1].set_xticklabels(statetypes_lab_main)
    for ax in axes[:-1]:
        ax.tick_params(axis='x', bottom=False, labelbottom=False)
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.02)
    plt.savefig((outpath+'/vert_SC_dFC_barplots.png'),dpi=720,bbox_inches='tight')
    plt.close()

    # Plot one plot for each cognitive variable. Plot bar plots for each R2,
    # where AC, MC, and S are grouped together. Modalities are spaced apart.
    # Annotate on top of the bar plots from top-down - the R2 value, then the percentages if 
    # they exist (on top of each other). No p-values.
    nstate = nstatetypes
    ncontrol = nctrl
    gap = 0.7
    barwidth = 0.8 / ncontrol
    xpos = []
    curr = 0
    for stateidx in range(nstate):
        for ctrlidx in range(ncontrol):
            xpos.append(curr)
            curr += 1
        curr += gap
    xpos = np.array(xpos)
    ctrl_colors = {'ave': '#1f77b4', 'mod': '#d62728', 'abs_deg': '#2ca02c'}
    fontsize_annot = 10
    fig, axes = plt.subplots(ncog, 1, figsize=(12,4.5*ncog),sharex=True,sharey=True)
    for cidx, ccog in enumerate(coglist):
        ax = axes[cidx]
        plotmn = plotmn_dict_supp[ccog]
        for ctrlidx, ctrltype in enumerate(controltypes):
            cvals = []
            barpos = []
            idx = 0
            for stateidx, statetype in enumerate(statetypes_supp):
                barpos.append(xpos[idx + ctrlidx])
                cvals.append(plotmn.loc[ctrltype, statetype])
                idx += ncontrol
            ax.bar(barpos, cvals, width=barwidth,
                    label=controltypes_lab[ctrlidx],
                    color=ctrl_colors[ctrltype],
                    edgecolor='black', linewidth=0.5)
            idx = 0
            for stateidx, statetype in enumerate(statetypes_supp):
                val = plotmn.loc[ctrltype, statetype]
                cx = xpos[idx + ctrlidx]
                idx += ncontrol
                bar_top = max(val, 0)
                x_offset = -9
                pt_offset = 3
                ckey = f'{ctrltype} {statetype} {ccog}'
                if statetype in statetype_perc_supp.keys():
                    percs = plotperc_dict_supp[ckey]
                    perclabs = statetype_perc_lab_supp[statetype]
                    for plabel, praw in zip(perclabs, percs.index):
                        ax.annotate(f'{plabel}: {percs[praw]:.1f}%', xy=(cx,bar_top), xycoords='data',
                                    xytext=(x_offset,pt_offset), textcoords='offset points',
                                    ha='left', va='bottom', fontsize=fontsize_annot)
                        pt_offset += 13
                ax.annotate(f'{val:.3f}', xy=(cx,bar_top), xycoords='data',
                            xytext=(x_offset,pt_offset), textcoords='offset points',
                            ha='left', va='bottom', fontsize=fontsize_annot, fontweight='bold')
        ax.axhline(0, color='black', linewidth=0.8)
        if cidx == 1:
            ax.set_ylabel('R²')
        ax.text(-0.10, 0.5, coglist_lab[cidx],transform=ax.transAxes,
            ha='center', va='center', rotation=0, fontsize=12, fontweight='bold')
        ax.set_ylim(*shared_ylim)
        if cidx == 0:
            ax.legend(title='Metric', loc='upper right')
    block_centers = []
    idx = 0
    for stateidx in range(nstate):
        block_xpos = xpos[idx:idx+ncontrol]
        block_centers.append(block_xpos.mean())
        idx += ncontrol
    axes[-1].set_xticks(block_centers)
    axes[-1].set_xticklabels(statetypes_lab_supp)
    for ax in axes[:-1]:
        ax.tick_params(axis='x', bottom=False, labelbottom=False)
    plt.tight_layout()
    fig.subplots_adjust(hspace=0.02)
    plt.savefig((outpath+'/vert_SC_sFC_dFC_barplots.png'),dpi=720,bbox_inches='tight')
    plt.close()
    print('Parameter variation done.')
