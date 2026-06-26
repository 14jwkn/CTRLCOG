# -*- coding: utf-8 -*-
"""
For a given k for clustering, type of SC normalization, type and percentage of thresholding, 
number of CV repetitions, inner k, outer k, analysis block, and number of permutations for the analysis, 
illustrate the relationship between Haufe scores and controllability or strength/principal cortical gradient using
scatter plots with the quadratic regression line, labelled with R2 and Spearman's correlation between Haufe scores 
and controllability or strength, with FDR-corrected significance color-coded.

Output:
ctrlpred_summary_plot.png Scatter plot with quadratic regression line and R2 and Spearman's correlation between Haufe scores and controllability or strength, with FDR-corrected significance color-coded.
grad1pred_summary_plot.png Scatter plot with quadratic regression line and R2 and Spearman's correlation between Haufe scores and the principal gradient, with FDR-corrected significance color-coded.

Usage: 
    LE_group_subcontrol_structfunc_KRR_maplook_plot.py <k> <sctype> <threstype> <thresval> <controltype> <statetype> <septype> <nrep> <inner_k> <outer_k> <wantblock> <wantperm>
    
Arguments:
    
    <k> k for k-means clustering
    <sctype> SC type of normalization
    <threstype> Threshold type
    <thresval> Threshold value
    <controltype> What control types are used
    <statetype> What statetypes are used
    <septype> What cognitive separation types are used
    <nrep> Number of CV repetitions
    <inner_k> Inner K in K-fold CV hyperparameter search
    <outer_k> Outer K in K-fold CV
    <wantblock> Which analysis block to run
    <wantperm> What number of permuations to run

"""

import os, sys, time, random
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from matplotlib import pyplot as plt
import matplotlib.lines as mlines
from docopt import docopt

#Quadratic regression.
def quad_reg(cx,cy):
   
    #Generate squared term and add to model.
    cx2 = cx**2
    full_x = np.column_stack([cx2,cx])

    #Fit model.
    cmodel = LinearRegression().fit(full_x,cy)

    #Extract coefficients and R2.
    y_pred = cmodel.predict(full_x)
    b2 = cmodel.coef_[0]
    b1 = cmodel.coef_[1]
    b0 = cmodel.intercept_
    r2 = cmodel.score(full_x,cy)

    #Return.
    return (b2,b1,b0,r2)

#Reader for control values.
def control_reader(nk,nroi,
                   basepath,scpath,sFCpath):

    #Read in the full Gramian error table to find states with errors.
    infile = (basepath+'/state_images/SC_sFC_dFC_gram.csv')
    gramall = pd.read_csv(infile,index_col=0,header=None)
    gramstates = gramall.index.values
    gramstates = gramstates[np.squeeze((gramall!=0).values)].tolist()
    ngram = len(gramstates)
        
    #Get FC AC, MC, and AD.
    infile = (basepath+'ave_tab.csv')
    dFCave = pd.read_csv(infile,index_col=0,header=None)
    infile = (basepath+'mod_tab.csv')
    dFCmod = pd.read_csv(infile,index_col=0,header=None)
    infile = (basepath+'abs_deg_tab.csv')
    dFCdeg = pd.read_csv(infile,index_col=0,header=None)

    #Generate FC predictor labels.
    predlabs = []
    for kidx in range(nk):
        for ridx in range(nroi):
            predlabs.append('s'+str(kidx+1)+'_r'+str(ridx+1))
    dFCave.columns = predlabs
    dFCmod.columns = predlabs
    dFCdeg.columns = predlabs

    #If the error list contains dFC states.
    for kidx in range(nk):
        errstate = ('s'+str(kidx+1))
        if (errstate in gramstates):
            
            #Read in the regions to remove.
            infile = (basepath+'/state_images/'+errstate+'_err.csv')
            errmat = pd.read_csv(infile,index_col=None,header=0)
            regrem = errmat.loc[:,'Region'].values.tolist()
            regrem = [(errstate+'_r'+str(x)) for x in regrem]
            newlabs = [x for x in predlabs if x not in regrem]

            #Remove those from the state.
            dFCave = dFCave.loc[:,newlabs]
            predlabs = dFCave.columns
            dFCmod = dFCmod.loc[:,newlabs]
            predlabs = dFCmod.columns
            dFCdeg = dFCdeg.loc[:,newlabs]
            predlabs = dFCdeg.columns

    #Read SC.
    infile = (scpath+'orig_ave_sc.csv')
    scave = pd.read_csv(infile,index_col=0,header=None)
    infile = (scpath+'orig_mod_sc.csv')
    scmod = pd.read_csv(infile,index_col=0,header=None)
    infile = (scpath+'deg_sc.csv')
    scdeg = pd.read_csv(infile,index_col=0,header=None)
    predlabs = [('sc_r'+str(ridx+1)) for ridx in range(nroi)]
    scave.columns = predlabs
    scmod.columns = predlabs
    scdeg.columns = predlabs

    #If the error list contains SC.
    errstate = 'SC'
    if (errstate in gramstates):
        
        #Read in the regions to remove.
        infile = (basepath+'/state_images/'+errstate+'_err.csv')
        errmat = pd.read_csv(infile,index_col=None,header=0)
        regrem = errmat.loc[:,'Region'].values.tolist()
        regrem = [('sc_r'+str(x)) for x in regrem]
        newlabs = [x for x in predlabs if x not in regrem]

        #Remove those from the state.
        scave = scave.loc[:,newlabs]
        scmod = scmod.loc[:,newlabs]
        scdeg = scdeg.loc[:,newlabs]

    #Read sFC.
    infile = (sFCpath+'ave_sFC.csv')
    sFCave = pd.read_csv(infile,index_col=0,header=None)
    infile = (sFCpath+'mod_sFC.csv')
    sFCmod = pd.read_csv(infile,index_col=0,header=None)
    infile = (sFCpath+'abs_deg_sFC.csv')
    sFCdeg = pd.read_csv(infile,index_col=0,header=None)
    predlabs = [('sFC_r'+str(ridx+1)) for ridx in range(nroi)]
    sFCave.columns = predlabs
    sFCmod.columns = predlabs
    sFCdeg.columns = predlabs

    #If the error list contains sFC.
    errstate = 'sFC'
    if (errstate in gramstates):
        
        #Read in the regions to remove.
        infile = (basepath+'/state_images/'+errstate+'_err.csv')
        errmat = pd.read_csv(infile,index_col=None,header=0)
        regrem = errmat.loc[:,'Region'].values.tolist()
        regrem = [('sFC_r'+str(x)) for x in regrem]
        newlabs = [x for x in predlabs if x not in regrem]

        #Remove those from the state.
        sFCave = sFCave.loc[:,newlabs]
        sFCmod = sFCmod.loc[:,newlabs]
        sFCdeg = sFCdeg.loc[:,newlabs]
        sFCgr = sFCgr.loc[newlabs]
    
    #Split single matrices into states.
    ctrl_collect = {}
    outmats = [scave,scmod,scdeg,
              sFCave,sFCmod,sFCdeg]
    outkeys = ['ave.sc','mod.sc','abs_deg.sc',
               'ave.sFC','mod.sFC','abs_deg.sFC']
    nout = len(outkeys)
    for ouidx in range(nout):
        outkey = outkeys[ouidx]
        outmat = outmats[ouidx]
        ctrl_collect[outkey] = outmat
   
    #Split multi-matrices into states.
    outmats = [dFCave,dFCmod,dFCdeg]
    outctrls = ['ave','mod','abs_deg']
    nout = len(outctrls)
    for cstate in [('s'+str(x+1)) for x in range(nk)]:
        for ouidx in range(nout):
            cctrl = outctrls[ouidx]
            inmat = outmats[ouidx]
            cpred = inmat.columns
            newpred = [x for x in cpred if cstate in x]
            outmat = inmat.loc[:,newpred]
            outkey = (cctrl+'.'+cstate)
            ctrl_collect[outkey] = outmat

    #Return each item.
    return (ctrl_collect)

#Reader for prediction values.
def predict_reader(basepath,nrep,inner_k,outer_k,sctype,
                   ctrltypes,nctrl,statetypes,nstatetype,septypes,nseptype,cog_septypes,
                   featver_list,nfeatver):
    
    #Produce predictive values for each model.
    pred_collect = {}
    for ctidx in range(nctrl):
        cctrl = ctrltypes[ctidx]

        #For each state type.
        for stidx in range(nstatetype):
            cstatetype = statetypes[stidx]

            #Set split state list based on state type.
            csplits = []
            if 'SC' in cstatetype:
                csplits = csplits + ['sc']
            if 'sFC' in cstatetype:
                csplits = csplits + ['sFC']
            if 'dFC' in cstatetype:
                csplits = csplits + [('s'+str(x+1)) for x in range(nk)]

            #For each cognitive set.
            for csidx in range(nseptype):
                cseptype = septypes[csidx]

                #Set cognitive list based on cognitive set.
                ccoglist = cog_septypes[csidx]
                cncog = len(ccoglist)

                #For each cognitive variable.
                for cidx in range(cncog):
                    ccog = ccoglist[cidx]
                
                    #For each feature type.
                    for fidx in range(nfeatver):
                        cfeatver = featver_list[fidx]

                        #Read in the version and average across folds and repetitions.
                        infile = (basepath+'KRRXFS/'+cctrl+'_'+cstatetype+'_'+cseptype+
                               nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+
                                '/score_collect.h5')
                        inkey = ('/rep_'+cfeatver+'_feat_'+ccog)
                        store = pd.HDFStore(infile,'r')
                        inmat = store.select(inkey)
                        store.close()
                        inmat = inmat.mean(axis=0)

                        #Split states.
                        for cstate in csplits:
                            cpred = inmat.index
                            newpred = [x for x in cpred if cstate in x]
                            outmat = inmat.loc[newpred]
                            outkey = (cctrl+'.'+cstatetype+'.'+ccog+'.'+cfeatver+'.'+cstate)
                            pred_collect[outkey] = outmat
    
    #Return each item.
    return pred_collect

if __name__ == '__main__':
    __spec__ = None
    
    #Catches arguments.
    args = docopt(__doc__)
    k = args['<k>']
    sctype = args['<sctype>']
    threstype = args['<threstype>']
    thresval = args['<thresval>']
    nrep = args['<nrep>']
    inner_k = args['<inner_k>']
    outer_k = args['<outer_k>']
    wantblock = args['<wantblock>']
    wantperm = args['<wantperm>']
    print('Doing:',k,sctype,threstype,thresval,
          nrep,inner_k,outer_k,wantblock,wantperm)
    
    #Set paths.
    controltype = 'controlcomp'
    statetype = 'statecomp'
    septype = 'cogcomp'
    subgroup = 'full'
    sc_subgroup = 'dr_full'
    basepath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                subgroup+'/'+k+'/SC_dFC/'+
                sc_subgroup+'/collect/'+threstype+'/'+thresval+'/')
    scpath = ('../outputs/d_SC/'+sc_subgroup+'/'+threstype+'/'+thresval+'/'+sctype+'/')
    sFCpath = ('../outputs/r_sFC/'+
                sc_subgroup+'/'+threstype+'/'+thresval+'/')
    explorepath = (basepath+'KRRXFS/'+controltype+'_'+statetype+'_'+septype+'_'+
                  nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/correxplore/TwoP_FDR/specialp/')
    os.makedirs(explorepath,exist_ok=True)
    
    #Initialize basic values.
    nroi = 360
    nk = int(k)
    xother_cont = ['Age']
    xother_cat_list = ['Gender']
    xother_cat = '|'.join(xother_cat_list)
    with open('dr_full_intersect.txt') as f:
        sublist = [label.rstrip() for label in f] 
    sublist = list(map(int,sublist))
    nsub = len(sublist)
    featver_list = ['covha']
    nfeatver = len(featver_list)
    predtype_list = ['full']
    npredtype = len(predtype_list)
    full_statesplits = ['sc','sFC'] + [('s'+str(x+1)) for x in range(nk)]

    #Set types of interest to be compared.
    singstates = ['SC']
    nsing = len(singstates)
    multstates = ['dFCcat','SC_dFCcat']
    nmult = len(multstates)
    statetypes = singstates + multstates
    nstatetype = len(statetypes)
    ctrltypes = ['ave','mod','abs_deg']
    nctrl = len(ctrltypes)
    septypes = ['comCFAng']
    nseptype = len(septypes)
    cog_septypes = [['gCFA','P24_CR','PV']]
    coglist = sum(cog_septypes,[])
    ncog = len(coglist)
    
    #Get cognitive and other predictor data.
    infile = ('../outputs/c_cognition/'+subgroup+'/pred_all.csv')
    othercog = pd.read_csv(infile,index_col=0,header=0)
    othercog = othercog.loc[sublist,:]
    cogmat = othercog.loc[:,coglist]
    othercont = othercog.loc[:,xother_cont]

    #Read in dummy data for categorical confounds, then merge.
    infile = ('../outputs/c_cognition/'+subgroup+'/dum_sel.csv')
    dummat = pd.read_csv(infile,index_col=0,header=0)
    dummat = dummat.loc[sublist,:]
    othercat = dummat.filter(regex=xother_cat)
    othermat = pd.concat((othercont,othercat),axis=1)

    #Read in control values.
    ctrl_collect = control_reader(nk,nroi,
                                  basepath,scpath,sFCpath)
    
    #Read in predictive values.
    pred_collect = predict_reader(basepath,nrep,inner_k,outer_k,sctype,
                                ctrltypes,nctrl,statetypes,nstatetype,septypes,nseptype,cog_septypes,
                                featver_list,nfeatver)
    
    #Read in gradient values.
    infile = ('../outputs/r_sFC/dr_full/none/0/sFC_gradients.csv')
    sFCgr_all = pd.read_csv(infile,header=None)
    infile = ('../outputs/r_sFC/dr_full/none/0/sFC_gradients_flip.csv')
    sFCgr_flip = pd.read_csv(infile,header=None).values.tolist()[0]
    ngr = len(sFCgr_flip)
    for gidx in range(ngr):
        if sFCgr_flip[gidx] == 'T':
            sFCgr_all.iloc[:,gidx] = -sFCgr_all.iloc[:,gidx]
    sFCgr_all.index = [('r'+str(ridx+1)) for ridx in range(nroi)]

    #Do group-average controllability relationship with predictiveness.
    print('Control-predictiveness relationship.')
    outlab = 'ctrlpred'
    outpath = (explorepath+outlab)

    #Set up.
    cstatetype = 'SC_dFCcat'
    csplits = []
    if 'SC' in cstatetype:
        csplits = csplits + ['sc']
    if 'sFC' in cstatetype:
        csplits = csplits + ['sFC']
    if 'dFC' in cstatetype:
        csplits = csplits + [('s'+str(x+1)) for x in range(nk)]
    nsplit = len(csplits)
    septypeint = 'comCFAng'
    cogint = ['PV','P24_CR','gCFA']
    ncogint = len(cogint)

    # Generate row-column combos.
    row_combos = []
    for ccog in cogint:
        if ccog == 'gCFA':
            c_ctrltypes = ['ave', 'mod', 'abs_deg']
        elif ccog == 'PV':
            c_ctrltypes = ['ave']
        elif ccog == 'P24_CR':
            c_ctrltypes = ['mod']
        for cctrl in c_ctrltypes:
            row_combos.append((ccog, cctrl))
    nrows = len(row_combos)     
    ncols = len(csplits)     

    # Gather scatter plots.
    fig, axes = plt.subplots(nrows,ncols,figsize=(3*ncols,3*nrows),squeeze=False)
    for row_idx,(ccog,cctrl) in enumerate(row_combos):

        # Convert to labels.
        if ccog == 'gCFA':
            ccog_lab = 'g'
        elif ccog == 'P24_CR':
            ccog_lab = 'gF'
        elif ccog == 'PV':
            ccog_lab = 'gC'
        if cctrl == 'ave':
            cctrl_lab = 'Average Controllability'
        elif cctrl == 'mod':
            cctrl_lab = 'Modal Controllability'
        elif cctrl == 'abs_deg':
            cctrl_lab = 'Strength'
            
        # Read quadratic regression and Spearman's correlation values.
        infile = (outpath+'/full/full_'+cctrl+'_'+cstatetype+
                '_comCFAng_covha_ctrlpred_'+ccog+'_summary.csv')
        inmat = pd.read_csv(infile,index_col=0)

        # Split into states.
        for col_idx,cstate in enumerate(csplits):

            # Select cell for plot.
            ax = axes[row_idx,col_idx]   

            # Get QR R2, QR R2 P, SP, SP P.
            cr2 = round(inmat.loc[cstate,'R2'],2)
            if (inmat.loc[cstate,'R2_P'] < 0.05):
                cr2_sig = '*'
            else:
                cr2_sig = ''
            csp = round(inmat.loc[cstate,'SP'],2)
            if (inmat.loc[cstate,'SP_P'] < 0.05):
                csp_sig = '*'
            else:
                csp_sig = ''
            
            # Turn title a color based on significance and sign.
            if (cr2_sig=='*') and (csp_sig=='*') and (csp >= 0):
                title_cl = '#FF0000'
            elif (cr2_sig=='*') and (csp_sig=='*') and (csp < 0):
                title_cl = '#00A2FF'
            elif (cr2_sig=='*') and (csp_sig!='*'):
                title_cl = '#097969'
            elif (cr2_sig!='*') and (csp_sig=='*'):
                title_cl = '#FFAA33'
            else:
                title_cl = 'black'
           
            # Get predictive score.
            inkey = (cctrl + '.' + cstatetype + '.' + ccog + '.covha.' + cstate)
            cpred = pred_collect[inkey]
            cpred_full = pd.Series(np.zeros(nroi),index=[(cstate+'_r'+str(x+1)) for x in range(nroi)])
            cpred_full[cpred.index] = cpred

            # Get controllability.
            inkey = (cctrl + '.' + cstate)
            cmat = ctrl_collect[inkey].mean()
            cmat_full = pd.Series(np.ones(nroi),index=[(cstate+'_r'+str(x+1)) for x in range(nroi)])
            cmat_full[cmat.index] = cmat

            # Concatenate and format.
            fullmat = pd.concat((cpred_full, cmat_full), axis=1)
            fullmat.columns = ['Predictive Score', 'Controllability']
            cX = fullmat.iloc[:,0]
            cY = fullmat.iloc[:,1]

            # Quadratic regression.
            beta = np.polyfit(cX, cY, deg=2)
            line_cX = np.linspace(min(cX), max(cX), 100)
            line_cY = np.polyval(beta, line_cX)

            # Plot scatter plot with values labelled.
            ax.scatter(cX,cY,s=2)
            ax.axvline(x=0, color='red', linestyle='--', linewidth=1)
            ax.plot(line_cX, line_cY, color='green')
            if cstate == 's3':
                xlab_out = f'{cctrl_lab} Regional Importance for Predicting {ccog_lab}'
            else:
                xlab_out = ''
            if cstate == 'sc':
                ylab_out = f'{cctrl_lab}'
            else:
                ylab_out = ''
            ax.set_xlabel(xlab_out)
            ax.set_ylabel(ylab_out)
            ax.set_title(
                f'{cstate.upper()} | QR R²={cr2}{cr2_sig} | RHO={csp}{csp_sig}',
                fontsize=10,loc='left',color=title_cl)

    # Draw a single line between the row above g, shift other plots down.
    plt.tight_layout()
    for i, (ccog, cctrl) in enumerate(row_combos):
        if ccog == 'gCFA':
            gCFA_start_row = i
            break
    gap = 0.03  
    for row_idx in range(gCFA_start_row, nrows):
        for col_idx in range(ncols):
            ax = axes[row_idx, col_idx]
            pos = ax.get_position()
            ax.set_position([pos.x0, pos.y0 - gap, pos.width, pos.height])
    ax_above = axes[gCFA_start_row - 1, 0]
    ax_below = axes[gCFA_start_row, 0]
    y_above = ax_above.get_position().y0
    y_below = ax_below.get_position().y1
    y_mid = (y_above + y_below) / 2
    line = mlines.Line2D(
        [0.01, 0.99], [y_mid, y_mid],
        transform=fig.transFigure,
        color='black', linewidth=1.5, linestyle='-'
    )
    fig.add_artist(line)

    # Save.
    plt.savefig((outpath+'/ctrlpred_summary_plot.png'),dpi=1080,bbox_inches='tight')
    plt.close()

    #Do sFC gradient 1 relationship with predictiveness.
    print('sFC gradient 1-predictiveness relationship.')
    outlab = 'grad1pred'
    outpath = (explorepath+outlab)

    #Set up.
    cstatetype = 'SC_dFCcat'
    csplits = []
    if 'SC' in cstatetype:
        csplits = csplits + ['sc']
    if 'sFC' in cstatetype:
        csplits = csplits + ['sFC']
    if 'dFC' in cstatetype:
        csplits = csplits + [('s'+str(x+1)) for x in range(nk)]
    nsplit = len(csplits)
    septypeint = 'comCFAng'
    cogint = ['gCFA']
    ncogint = len(cogint)
    sFCgr = sFCgr_all.iloc[:,0]

    # Generate row-column combos.
    row_combos = []
    if ccog == 'gCFA':
        c_ctrltypes = ['ave', 'mod', 'abs_deg']
    for cctrl in c_ctrltypes:
        row_combos.append((ccog, cctrl))
    nrows = len(row_combos)     
    ncols = len(csplits)     

    # Gather scatter plots.
    fig, axes = plt.subplots(nrows,ncols,figsize=(3*ncols,3*nrows),squeeze=False)
    for row_idx,(ccog,cctrl) in enumerate(row_combos):

        # Convert to labels.
        if ccog == 'gCFA':
            ccog_lab = 'g'
        if cctrl == 'ave':
            cctrl_lab = 'Average Controllability'
        elif cctrl == 'mod':
            cctrl_lab = 'Modal Controllability'
        elif cctrl == 'abs_deg':
            cctrl_lab = 'Strength'
            
        # Read quadratic regression and Spearman's correlation values.
        infile = (outpath+'/full/full_'+cctrl+'_'+cstatetype+
                '_comCFAng_covha_grad1pred_'+ccog+'_summary.csv')
        inmat = pd.read_csv(infile,index_col=0)

        # Split into states.
        for col_idx,cstate in enumerate(csplits):

            # Select cell for plot.
            ax = axes[row_idx,col_idx]   

            # Get QR R2, QR R2 P, SP, SP P.
            cr2 = round(inmat.loc[cstate,'R2'],2)
            if (inmat.loc[cstate,'R2_P'] < 0.05):
                cr2_sig = '*'
            else:
                cr2_sig = ''
            csp = round(inmat.loc[cstate,'SP'],2)
            if (inmat.loc[cstate,'SP_P'] < 0.05):
                csp_sig = '*'
            else:
                csp_sig = ''
            
            # Turn title a color based on significance and sign.
            if (cr2_sig=='*') and (csp_sig=='*') and (csp >= 0):
                title_cl = '#FF0000'
            elif (cr2_sig=='*') and (csp_sig=='*') and (csp < 0):
                title_cl = '#00A2FF'
            elif (cr2_sig=='*') and (csp_sig!='*'):
                title_cl = '#097969'
            elif (cr2_sig!='*') and (csp_sig=='*'):
                title_cl = "#FFAA33"
            else:
                title_cl = 'black'
            
            # Get regional importances.
            inkey = (cctrl+'.'+cstatetype+'.'+ccog+'.covha.'+cstate)
            cpred = pred_collect[inkey]
            
            #Requires full 360 regions. For predictiveness, replace with 0 for zero predictiveness.
            cpred_full = pd.Series(np.zeros(nroi),index=[(cstate+'_r'+str(x+1)) for x in range(nroi)])
            cpred_full[cpred.index] = cpred

            #Index comparator. Get all original values.
            cmat = sFCgr.copy()
            cmat.index = [(cstate+'_'+x) for x in cmat.index]
            cmat_full = cmat
            
            #Combine.
            fullmat = pd.concat((cpred_full,cmat_full),axis=1)
            fullmat.columns = ['Predictive Score','Gradient 1']
            cX = fullmat.iloc[:,0]
            cY = fullmat.iloc[:,1]

            # Quadratic regression.
            beta = np.polyfit(cX, cY, deg=2)
            line_cX = np.linspace(min(cX), max(cX), 100)
            line_cY = np.polyval(beta, line_cX)

            # Plot scatter plot with values labelled.
            ax.scatter(cX,cY,s=2)
            ax.axvline(x=0, color='red', linestyle='--', linewidth=1)
            ax.plot(line_cX, line_cY, color='green')
            if cstate == 's3':
                xlab_out = f'{cctrl_lab} Regional Importance for Predicting {ccog_lab}'
            else:
                xlab_out = ''
            if cstate == 'sc':
                ylab_out = f'Principal Gradient'
            else:
                ylab_out = ''
            ax.set_xlabel(xlab_out)
            ax.set_ylabel(ylab_out)
            ax.set_title(
                f'{cstate.upper()} | QR R²={cr2}{cr2_sig} | RHO={csp}{csp_sig}',
                fontsize=10,loc='left',color=title_cl)

    # Shift plots down.
    plt.tight_layout()
    for i, (ccog, cctrl) in enumerate(row_combos):
        if ccog == 'gCFA':
            gCFA_start_row = i
            break
    gap = 0.03  
    for row_idx in range(gCFA_start_row, nrows):
        for col_idx in range(ncols):
            ax = axes[row_idx, col_idx]
            pos = ax.get_position()
            ax.set_position([pos.x0, pos.y0 - gap, pos.width, pos.height])
            
    # Save.
    plt.savefig((outpath+'/grad1pred_summary_plot.png'),dpi=1080,bbox_inches='tight')
    plt.close()
    print('Exploration done.')
