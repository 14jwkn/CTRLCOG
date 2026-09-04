# -*- coding: utf-8 -*-
"""
For a given number of repetitions, inner loop number, and outer loop number for CV, 
for the main model, plot the percentage contributions of each dFC state separately
along with SC in for each cognitive variable.
Output:
all_percsep.png Line plots for each cognitive variable for AC, MC, and S for matrix-wise percentage contributions.

Usage: 
    LE_group_subcontrol_structfunc_KRR_score_percsep.py <nrep> <inner_k> <outer_k>
    
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
    statetype = 'SC_dFCsep'
    threstype = 'groupconsist'
    thresval = '50'
    septype = 'comCFAng'
    
    # Plot one plot for each cognitive variable, side by side. Each plot contains
    # line plots where AC, MC, and S are different colors. All plots have the same axes.
    nk = int(k)
    statelist = ['sc'] + [f's{x+1}' for x in range(nk)]
    statelist_lab = ['SC'] + [f'S{x+1}' for x in range(nk)]
    nstate = len(statelist)
    controltypes = ['ave','mod','abs_deg'] 
    controltypes_lab = ['AC','MC','S']
    nctrl = len(controltypes)
    coglist = ['gCFA','P24_CR','PV']
    coglist_lab = ['g','gF','gC']
    ncog = len(coglist)
    plotmat_list = []
    for ccog in coglist:
        plotmat = pd.DataFrame(np.zeros((nctrl,nstate)),
                               index=controltypes,columns=statelist)
        for ctrltype in controltypes:

            # Read.
            inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'
                    + subgroup + '/'
                    + k + '/SC_dFC/' + sc_subgroup + '/collect/'
                    + threstype + '/' + thresval + '/KRRXFS/'
                    + ctrltype + '_' + statetype + '_' + septype + '_'
                    + nrep + '_' + inner_k + '_'
                    + outer_k + '_' + sctype)
            infile = (inpath + '/splitacc_strict_perc.csv')
            plotmat.loc[ctrltype,:] = pd.read_csv(infile,index_col=0).loc[statelist,ccog]

        # Append.
        plotmat_list.append(plotmat)
        print(plotmat.to_string())

    # Plot.
    ctrl_colors = {'ave': '#1f77b4', 'mod': '#d62728', 'abs_deg': '#2ca02c'}
    xpos = np.arange(nstate)
    fig, axes = plt.subplots(1, ncog, figsize=(6*ncog, 5), sharey=True, sharex=True)
    for cidx, ccog in enumerate(coglist):
        ax = axes[cidx]
        plotmat = plotmat_list[cidx]
        for ctrltype in controltypes:
            ax.plot(xpos, plotmat.loc[ctrltype, :].values,
                    label=controltypes_lab[controltypes.index(ctrltype)],
                    color=ctrl_colors[ctrltype])
        ax.set_xticks(xpos)
        ax.set_xticklabels(statelist_lab)
        ax.set_title(coglist_lab[cidx])
        if cidx == 1:
            ax.set_xlabel('State')
        if cidx == 0:
            ax.set_ylabel('R² Percentage')
        else:
            ax.tick_params(axis='y', which='both', left=False)
    handles, labels = axes[-1].get_legend_handles_labels()
    fig.legend(handles, labels, title='Metric', loc='upper left',
               bbox_to_anchor=(1.0, 0.94))
    plt.tight_layout()
    outfile = (outpath+'/all_percsep.png')
    plt.savefig(outfile, dpi=720, bbox_inches='tight')
    plt.close()
    print('Saved.')
