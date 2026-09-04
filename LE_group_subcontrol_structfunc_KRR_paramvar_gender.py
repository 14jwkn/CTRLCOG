# -*- coding: utf-8 -*-
"""
For a given number of repetitions, inner loop number, and outer loop number for CV, 
find out whether the main findings of the paper hold across different parameter variations
for the SC & dFC model - dFC performs better than SC and controllability and strength exhibit
high correspondence, g is predicted better than gF and gC, and the principal gradient exhibits a relationship 
with regional importance in the expected direction for dFC. The parameters that is varied is gender. Generate line plots 
that display these comparisons.
Output:
gender_SC_vs_dFC.png Compare dFC and SC, and controllability and strength, across genders.
gender_g_vs_gF_vs_gC.png Compare g and gF and gC across genders.
gender_sFCgr.png  Compare the principal gradient and regional importance across genders.

Usage: 
    LE_group_subcontrol_structfunc_KRR_paramvar_gender.py <nrep> <inner_k> <outer_k>
    
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

    # Set parameters.
    k = '6'
    sctype = 'fpt'
    threstype = 'groupconsist'
    thresval = '50'
    septype = 'comCFAng'
    nrep = '5'
    inner_k = '5'
    outer_k = '10'
    genders = ['msub','fsub']
    genders_lab = ['Male','Female']
    
    # --------------------------------- SC vs dFC -------------------------------- #

    # Set parameters.
    inkey = '/rep_testacc'
    controltypes = ['ave','mod','abs_deg'] 
    statetypes = ['SC_dFCcat','dFCcat','SC']
    controltypes_lab = ['AC','MC','S']
    statetypes_lab = ['SC & dFC','dFC','SC']
    ccog = 'gCFA'

    # Plot the prediction of g by each as bar plots, male and female separated.
    # Each state type is also separated by a gap.
    complist = [f'{x}_{y}' for x in statetypes for y in controltypes]
    ncomp = len(complist)
    plotmat = pd.DataFrame(np.zeros((ncomp,2)),index=complist,columns=genders)
    for gender in genders:
        for statetype in statetypes:
            for ctrltype in controltypes:
                clab = f'{statetype}_{ctrltype}'
                inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                        + subgroup + '/'
                        + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                        + threstype + '/' + thresval + '/KRRXFS/'+gender+'_'
                        + ctrltype + '_' + statetype + '_' + septype + '_'
                        + nrep + '_' + inner_k + '_'
                        + outer_k + '_' + sctype)
                infile = (inpath + '/score_collect.h5')
                store = pd.HDFStore(infile,'r')
                testacc = store.select(inkey)
                store.close()
                plotmat.loc[clab,gender] = testacc.mean(axis=0).loc[ccog]

    # Plot.
    ncontrol = len(controltypes)
    nstate = len(statetypes)
    gap = 0.7  
    barwidth = 0.35
    xpos = []
    curr = 0
    for stateidx in range(nstate):
        for ctrlidx in range(ncontrol):
            xpos.append(curr)
            curr += 1
        curr += gap  
    xpos = np.array(xpos)
    fig, ax = plt.subplots(figsize=(8,6))
    colors = {'msub': '#4C72B0','fsub': '#DD8452'}
    for gidx, gender in enumerate(genders):
        offset = (gidx - 0.5) * barwidth
        ax.bar(xpos + offset, plotmat[gender].values, width=barwidth,
               label=genders_lab[gidx], color=colors[gender])
    ax.set_xticks(xpos)
    ax.set_xticklabels(controltypes_lab*nstate)
    block_centers = []
    idx = 0
    for stateidx in range(nstate):
        block_xpos = xpos[idx:idx+ncontrol]
        block_centers.append(block_xpos.mean())
        idx += ncontrol
    for center, statetype_lab in zip(block_centers, statetypes_lab):
        ax.text(center, ax.get_ylim()[0] - 0.08*(ax.get_ylim()[1]-ax.get_ylim()[0]),
                statetype_lab, ha='center', va='top', fontsize=11)
    ax.set_ylabel('R²')
    ax.legend(title='Gender')
    ax.axhline(0, color='black', linewidth=0.8)
    plt.tight_layout()
    plt.savefig((outpath+'/gender_SC_vs_dFC.png'),dpi=720,bbox_inches='tight')
    plt.close()

    # ------------------------------- g vs gF vs gC ------------------------------ #

    # Set parameters.
    inkey = '/rep_testacc'
    controltypes = ['ave','mod','abs_deg'] 
    controltypes_lab = ['AC','MC','S']
    statetype = 'SC_dFCcat'
    statetype_lab = 'SC & dFC'
    coglist = ['gCFA','P24_CR','PV']
    coglist_lab = ['g','gF','gC']

    # Plot the prediction of g, gF, and gC by each as bar plots, male and female separated.
    # Each cognitive variable is also separated by a gap.
    complist = [f'{x}_{y}' for x in coglist for y in controltypes]
    ncomp = len(complist)
    plotmat = pd.DataFrame(np.zeros((ncomp,2)),index=complist,columns=genders)
    for gender in genders:
        for ccog in coglist:
            for ctrltype in controltypes:
                clab = f'{ccog}_{ctrltype}'
                inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                        + subgroup + '/'
                        + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                        + threstype + '/' + thresval + '/KRRXFS/'+gender+'_'
                        + ctrltype + '_' + statetype + '_' + septype + '_'
                        + nrep + '_' + inner_k + '_'
                        + outer_k + '_' + sctype)
                infile = (inpath + '/score_collect.h5')
                store = pd.HDFStore(infile,'r')
                testacc = store.select(inkey)
                store.close()
                plotmat.loc[clab,gender] = testacc.mean(axis=0).loc[ccog]

    # Plot.
    ncontrol = len(controltypes)
    ncog = len(coglist)
    gap = 0.7
    barwidth = 0.35
    xpos = []
    curr = 0
    for cogidx in range(ncog):
        for ctrlidx in range(ncontrol):
            xpos.append(curr)
            curr += 1
        curr += gap  
    xpos = np.array(xpos)
    fig, ax = plt.subplots(figsize=(8,6))
    colors = {'msub': '#4C72B0','fsub': '#DD8452'}
    for gidx, gender in enumerate(genders):
        offset = (gidx - 0.5) * barwidth
        ax.bar(xpos + offset, plotmat[gender].values, width=barwidth,
               label=genders_lab[gidx], color=colors[gender])
    ax.set_xticks(xpos)
    ax.set_xticklabels(controltypes_lab*ncog)
    block_centers = []
    idx = 0
    for cogidx in range(ncog):
        block_xpos = xpos[idx:idx+ncontrol]
        block_centers.append(block_xpos.mean())
        idx += ncontrol
    for center, cog_lab in zip(block_centers, coglist_lab):
        ax.text(center, ax.get_ylim()[0] - 0.08*(ax.get_ylim()[1]-ax.get_ylim()[0]),
                cog_lab, ha='center', va='top', fontsize=11)
    ax.set_ylabel('R²')
    ax.legend(title='Gender')
    ax.axhline(0, color='black', linewidth=0.8)
    plt.tight_layout()
    plt.savefig((outpath+'/gender_g_vs_gF_vs_gC.png'),dpi=720,bbox_inches='tight')
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

    # Set parameters.
    inkey = ('/rep_covha_feat_gCFA')
    controltypes = ['ave','mod','abs_deg'] 
    statetype = 'SC_dFCcat'
    controltypes_lab = ['AC','MC','S']
    statetype_lab = 'SC & dFC'
    ccog = 'gCFA'
    nk = int(k)
    klabs = [f's{str(x+1)}' for x in range(nk)] 
    klabs_lab = [f'S{str(x+1)}' for x in range(nk)] 

    # Get each dFC state and find correlation with gradient. Plot bar plots for 
    # the genders separately.
    complist = controltypes
    ncomp = len(complist)
    mmat = pd.DataFrame(np.zeros((ncomp,nk)),index=complist,columns=klabs)
    fmat = pd.DataFrame(np.zeros((ncomp,nk)),index=complist,columns=klabs)
    for ctrltype in controltypes:
        clab = ctrltype
        for gender in genders:
            inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                    + subgroup + '/'
                    + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                    + threstype + '/' + thresval + '/KRRXFS/'+gender+'_'
                    + ctrltype + '_' + statetype + '_' + septype + '_'
                    + nrep + '_' + inner_k + '_'
                    + outer_k + '_' + sctype)
            infile = (inpath + '/score_collect.h5')
            store = pd.HDFStore(infile,'r')
            featimp = store.select(inkey)
            store.close()
            allfeat = featimp.mean(axis=0)
            for kidx in range(nk):
                cstate = klabs[kidx]
                fcols = [x for x in allfeat.index if f'{cstate}_' in x]
                cstate_feat = allfeat[fcols]
                gcols = [x.replace(f'{cstate}_','') for x in fcols]
                cgr = sFCgr[gcols]
                if gender == 'msub':
                    mmat.loc[clab,cstate] = spearmanr(cstate_feat,cgr).statistic
                elif gender == 'fsub':
                    fmat.loc[clab,cstate] = spearmanr(cstate_feat,cgr).statistic

    # Plot.
    ncontrol = len(controltypes)
    gap = 0.7
    barwidth = 2 / ncontrol
    xpos = []
    curr = 0
    for kidx in range(nk):
        for ctrlidx in range(ncontrol):
            xpos.append(curr)
            curr += 1
        curr += gap  
    xpos = np.array(xpos)
    ctrl_colors = {'ave': '#1f77b4', 'mod': '#d62728', 'abs_deg': '#2ca02c'}
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    mats = {'msub': mmat, 'fsub': fmat}
    for aidx, gender in enumerate(genders):
        ax = axes[aidx]
        cmat = mats[gender]
        for cidx, ctrltype in enumerate(controltypes):
            cvals = []
            barpos = []
            idx = 0
            for kidx in range(nk):
                barpos.append(xpos[idx + cidx])
                cvals.append(cmat.loc[ctrltype, klabs[kidx]])
                idx += ncontrol
            ax.bar(barpos, cvals, width=barwidth, label=controltypes_lab[cidx],
                   color=ctrl_colors[ctrltype])
        block_centers = []
        idx = 0
        for kidx in range(nk):
            block_xpos = xpos[idx:idx+ncontrol]
            block_centers.append(block_xpos.mean())
            idx += ncontrol
        ax.set_xticks(block_centers)
        ax.set_xticklabels(klabs_lab, fontsize=11)
        ax.axhline(0, color='black', linewidth=0.8)
        ax.set_title(genders_lab[aidx], fontsize=13)
        ax.set_xlabel('State', fontsize=11)
        if aidx == 0:
            ax.set_ylabel('Principal Gradient RHO', fontsize=11)
        ax.tick_params(axis='y', labelsize=10)
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, title='Metrics', loc='upper left',
               bbox_to_anchor=(1.0, 1.0), fontsize=10, title_fontsize=11)
    plt.tight_layout()
    plt.savefig((outpath+'/gender_sFCgr.png'),dpi=720,bbox_inches='tight')
    plt.close()
    print('Parameter variation done.')
