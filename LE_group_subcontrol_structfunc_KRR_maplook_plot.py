# -*- coding: utf-8 -*-
"""
For a given k for clustering, type of SC normalization, and type and percentage of thresholding, 
illustrate the relationship between Haufe scores and controllability/principal cortical gradient using
linear and quadratic regression lines and Pearson's and Spearman's correlation.
Output:
ccog+'_'+cctrl+'_'+cstate+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'.jpg' Scatterplot with linear regression line and Pearson's correlation labelled.
quad_'+ccog+'_'+cctrl+'_'+cstate+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_covha_'+outlab+'.jpg' Scatterplot with quadratic regression line and Spearman's correlation labelled.
cstatetype+'_'+septypeint+'_'+corrtype+'_covha_corr.csv' Pearson's or Spearman's correlations collected.

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
from docopt import docopt

#Reader for control and degree values.
def control_reader(nk,nroi,basepath,scpath,sFCpath):

    #Read in the full Gramian error table to find states with errors.
    infile = (basepath+'/state_images/SC_sFC_dFC_gram.csv')
    gramall = pd.read_csv(infile,index_col=0,header=None)
    gramstates = gramall.index.values
    gramstates = gramstates[np.squeeze((gramall!=0).values)].tolist()
    ngram = len(gramstates)
        
    #Get FC AC, MC, and D.
    infile = (basepath+'ave_tab.csv')
    dFCave = pd.read_csv(infile,index_col=0,header=None)
    infile = (basepath+'mod_tab.csv')
    dFCmod = pd.read_csv(infile,index_col=0,header=None)
    infile = (basepath+'deg_tab.csv')
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
    infile = (sFCpath+'deg_sFC.csv')
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
    outkeys = ['ave.sc','mod.sc','deg.sc',
               'ave.sFC','mod.sFC','deg.sFC']
    nout = len(outkeys)
    for ouidx in range(nout):
        outkey = outkeys[ouidx]
        outmat = outmats[ouidx]
        ctrl_collect[outkey] = outmat
   
    #Split multi-matrices into states.
    outmats = [dFCave,dFCmod,dFCdeg]
    outctrls = ['ave','mod','deg']
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

#Reader for prediction accuracy values and p-values.
def predict_reader(basepath,nrep,inner_k,outer_k,
                   sctype,ctrltypes,nctrl,statetypes,nstatetype,septypes,nseptype,cog_septypes,
                   pstatetypes,npstatetype,pseptypes,npseptype,pcog_septypes):
    
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
                
                    #Read in the version and average across folds and repetitions.
                    infile = (basepath+'KRRXFS/'+cctrl+'_'+cstatetype+'_'+cseptype+
                                '_'+nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+
                                '/score_collect.h5')
                    inkey = ('/rep_covha_feat_'+ccog)
                    store = pd.HDFStore(infile,'r')
                    inmat = store.select(inkey)
                    store.close()
                    inmat = inmat.mean(axis=0)

                    #Split states.
                    for cstate in csplits:
                        cpred = inmat.index
                        newpred = [x for x in cpred if cstate in x]
                        outmat = inmat.loc[newpred]
                        outkey = (cctrl+'.'+cstatetype+'.'+ccog+'.'+cstate)
                        pred_collect[outkey] = outmat

    #Produce predictive values for each model.
    pval_collect = {}
    for ctidx in range(nctrl):
        cctrl = ctrltypes[ctidx]

        #For each state type.
        for stidx in range(npstatetype):
            cstatetype = pstatetypes[stidx]

            #Set split state list based on state type.
            csplits = []
            if 'SC' in cstatetype:
                csplits = csplits + ['sc']
            if 'sFC' in cstatetype:
                csplits = csplits + ['sFC']
            if 'dFC' in cstatetype:
                csplits = csplits + [('s'+str(x+1)) for x in range(nk)]

            #For each cognitive set.
            for csidx in range(npseptype):
                cseptype = pseptypes[csidx]

                #Set cognitive list based on cognitive set.
                ccoglist = pcog_septypes[csidx]
                cncog = len(ccoglist)

                #For each cognitive variable.
                for cidx in range(cncog):
                    ccog = ccoglist[cidx]
                
                    #Read in the version.
                    inpath = (basepath+'KRRXFS/controlcomp_statecomp_cogcomp'+
                                '_'+nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/')
                    infile = (inpath+cctrl+'_'+cstatetype+'_'+cseptype+'_TwoP_BHFDR_covha_feat.csv')
                    inmat = pd.read_csv(infile,index_col=0)
                    inmat = inmat.loc[:,ccog]
                    
                    #Split states.
                    for cstate in csplits:
                        cpred = inmat.index
                        newpred = [x for x in cpred if cstate in x]
                        outmat = inmat.loc[newpred]
                        outkey = (cctrl+'.'+cstatetype+'.'+ccog+'.'+cstate)
                        pval_collect[outkey] = outmat < 0.05
    
    #Return each item.
    return (pred_collect,pval_collect)

if __name__ == '__main__':
    __spec__ = None
    
    #Catches arguments.
    args = docopt(__doc__)
    k = args['<k>']
    sctype = args['<sctype>']
    threstype = args['<threstype>']
    thresval = args['<thresval>']
    controltype = args['<controltype>']
    statetype = args['<statetype>']
    septype = args['<septype>']
    nrep = args['<nrep>']
    inner_k = args['<inner_k>']
    outer_k = args['<outer_k>']
    wantblock = args['<wantblock>']
    wantperm = args['<wantperm>']
    print('Doing:',k,sctype,threstype,thresval,
          controltype,statetype,septype,
          nrep,inner_k,outer_k,wantblock,wantperm)
    
    #Set paths.
    subgroup = 'full'
    sc_subgroup = 'dr_full'
    basepath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                subgroup+'/'+k+'/SC_dFC/'+
                sc_subgroup+'/collect/'+threstype+'/'+thresval+'/')
    scpath = ('../outputs/d_SC/'+sc_subgroup+'/'+threstype+'/'+thresval+'/'+sctype+'/')
    sFCpath = ('../outputs/r_sFC/'+
                sc_subgroup+'/'+threstype+'/'+thresval+'/')
    explorepath = (basepath+'KRRXFS/'+controltype+'_'+statetype+'_'+septype+'_'+
                  nrep+'_'+inner_k+
                  '_'+outer_k+'_'+sctype+'/correxplore/plot/')
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
    full_statesplits = ['sc','sFC'] + [('s'+str(x+1)) for x in range(nk)]

    #Set types of interest to be compared.
    singstates = ['SC']
    nsing = len(singstates)
    multstates = ['dFCcat','SC_dFCcat']
    nmult = len(multstates)
    statetypes = singstates + multstates
    nstatetype = len(statetypes)
    ctrltypes = ['ave','mod']
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
    (pred_collect,pval_collect) = predict_reader(basepath,nrep,inner_k,outer_k,sctype,
                                                ctrltypes,nctrl,statetypes,nstatetype,septypes,nseptype,cog_septypes)
    
    #Read in gradient values.
    infile = ('../outputs/r_sFC/unclean/both/whole/dr_full/none/0/sFC_gradients.csv')
    sFCgr_all = pd.read_csv(infile,header=None)
    infile = ('../outputs/r_sFC/unclean/both/whole/dr_full/none/0/sFC_gradients_flip.csv')
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
    os.makedirs(outpath,exist_ok=True)

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
    cogint = ['gCFA','P24_CR','PV']
    ncogint = len(cogint)

    #Set up to save correlations.
    pr_corrsave_all = []
    sp_corrsave_all = []
    corrsave_labs = []
    for ccog in cogint:

        #Set control types of interest for the cognitive variable.
        if ccog=='gCFA':
            c_ctrltypes = ['ave','mod']
        elif ccog=='P24_CR':
            c_ctrltypes = ['mod']
        elif ccog=='PV':
            c_ctrltypes = ['ave']

        #Do comparisons.
        for cctrl in c_ctrltypes:

            #Set up to save correlations.
            corrsave_labs.append(ccog+' '+cctrl)
            pr_corrsave = []
            sp_corrsave = []

            #Do each state.
            for cstate in csplits:

                #Index predictiveness.
                inkey = (cctrl+'.'+cstatetype+'.'+ccog+'.'+cstate)
                cpred = pred_collect[inkey]
                
                #Requires full 360 regions. For predictiveness, replace with 0 for zero predictiveness.
                cpred_full = pd.Series(np.zeros(nroi),index=[(cstate+'_r'+str(x+1)) for x in range(nroi)])
                cpred_full[cpred.index] = cpred
                
                #Index comparator, generate group-average.
                inkey = (cctrl+'.'+cstate)
                cmat = ctrl_collect[inkey].mean()

                #Requires full 360 regions. For control, replace with 1 for what the error value was.
                cmat_full = pd.Series(np.ones(nroi),index=[(cstate+'_r'+str(x+1)) for x in range(nroi)])
                cmat_full[cmat.index] = cmat
                
                #Combine.
                fullmat = pd.concat((cpred_full,cmat_full),axis=1)
                fullmat.columns = ['Predictive Score','Controllability']
                cX = fullmat.iloc[:,0]
                cY = fullmat.iloc[:,1]

                #Generate Pearson's correlation.
                pr_corr = pearsonr(cX,cY).statistic
                corr_r = round(pr_corr,2)
                if corr_r == 0:
                    corr_r = 0
                
                #Add.
                pr_corrsave.append(corr_r)
                
                #Plot scatterplot.
                plt.scatter(cX,cY,s=2)
                plt.axvline(x=0,color='red',linestyle='--', linewidth=1)
                plt.xlabel(fullmat.columns[0])
                plt.ylabel(fullmat.columns[1])
                plt.title((ccog+' '+cctrl+' '+cstate+' | r = '+str(corr_r)),loc='right')

                #Save.
                corrtype = 'pearson'
                outfile = (outpath+'/'+ccog+'_'+cctrl+'_'+cstate+'_'+
                        cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'.jpg')
                plt.savefig(outfile,bbox_inches='tight',dpi=360)  
                plt.close()

                #Do quadratic regression.
                beta = np.polyfit(cX,cY,deg=2)
                line_cX = np.linspace(min(cX),max(cX),100)
                line_cY = np.polyval(beta,line_cX)

                #Generate Spearman's correlation.
                sp_corr = spearmanr(cX,cY).statistic
                corr_r = round(sp_corr,2)
                if corr_r == 0:
                    corr_r = 0
                
                #Add.
                sp_corrsave.append(corr_r)
                
                #Plot scatterplot.
                plt.scatter(cX,cY,s=2)
                plt.axvline(x=0,color='red',linestyle='--', linewidth=1)
                plt.plot(line_cX,line_cY,color='green')
                plt.xlabel(fullmat.columns[0])
                plt.ylabel(fullmat.columns[1])
                plt.title((ccog+' '+cctrl+' '+cstate+' | SP r = '+str(corr_r)),loc='right')

                #Save.
                corrtype = 'spearman'
                outfile = (outpath+'/quad_'+ccog+'_'+cctrl+'_'+cstate+'_'+
                        cstatetype+'_'+septypeint+'_'+corrtype+'_covha_'+outlab+'.jpg')
                plt.savefig(outfile,bbox_inches='tight',dpi=360)  
                plt.close()

            #Add.
            pr_corrsave_all.append(pr_corrsave)
            sp_corrsave_all.append(sp_corrsave)
    
    #Format and save.
    corrtype = 'pearson'
    corrsave_out = pd.DataFrame(pr_corrsave_all,index=corrsave_labs,columns=csplits)
    outfile = (outpath+'/'+cstatetype+'_'+septypeint+'_'+corrtype+'_covha_corr.csv')
    corrsave_out.to_csv(outfile,index=True,header=True)
    corrtype = 'spearman'
    corrsave_out = pd.DataFrame(sp_corrsave_all,index=corrsave_labs,columns=csplits)
    outfile = (outpath+'/'+cstatetype+'_'+septypeint+'_'+corrtype+'_covha_corr.csv')
    corrsave_out.to_csv(outfile,index=True,header=True)
    
    #Do sFC gradient 1 relationship with predictiveness.
    print('sFC gradient 1-predictiveness relationship.')
    outlab = 'grad1pred'
    outpath = (explorepath+outlab)
    os.makedirs(outpath,exist_ok=True)

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

    #Set up to save correlations.
    pr_corrsave_all = []
    sp_corrsave_all = []
    corrsave_labs = []
    for ccog in cogint:

        #Set control types of interest for the cognitive variable.
        if ccog=='gCFA':
            c_ctrltypes = ['ave','mod']
        elif ccog=='P24_CR':
            c_ctrltypes = ['mod']
        elif ccog=='PV':
            c_ctrltypes = ['ave']

        #Do comparisons.
        for cctrl in c_ctrltypes:

            #Set up to save correlations.
            corrsave_labs.append(ccog+' '+cctrl)
            pr_corrsave = []
            sp_corrsave = []

            #Do each state.
            for cstate in csplits:

                #Index predictiveness.
                inkey = (cctrl+'.'+cstatetype+'.'+ccog+'.'+cstate)
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

                #Generate Pearson's correlation.
                pr_corr = pearsonr(cX,cY).statistic
                corr_r = round(pr_corr,2)
                if corr_r == 0:
                    corr_r = 0
                
                #Add.
                pr_corrsave.append(corr_r)
                
                #Plot scatterplot.
                plt.scatter(cX,cY,s=2)
                plt.axvline(x=0,color='red',linestyle='--', linewidth=1)
                plt.xlabel(fullmat.columns[0])
                plt.ylabel(fullmat.columns[1])
                plt.title((ccog+' '+cctrl+' '+cstate+' | r = '+str(corr_r)),loc='right')

                #Save.
                corrtype = 'pearson'
                outfile = (outpath+'/'+ccog+'_'+cctrl+'_'+cstate+'_'+
                        cstatetype+'_'+septypeint+'_'+corrtype+'_covha_'+outlab+'.jpg')
                plt.savefig(outfile,bbox_inches='tight',dpi=360)  
                plt.close()

                #Do quadratic regression.
                beta = np.polyfit(cX,cY,deg=2)
                line_cX = np.linspace(min(cX),max(cX),100)
                line_cY = np.polyval(beta,line_cX)

                #Generate Spearman's correlation.
                sp_corr = spearmanr(cX,cY).statistic
                corr_r = round(sp_corr,2)
                if corr_r == 0:
                    corr_r = 0
                
                #Add.
                sp_corrsave.append(corr_r)
                
                #Plot scatterplot.
                plt.scatter(cX,cY,s=2)
                plt.axvline(x=0,color='red',linestyle='--', linewidth=1)
                plt.plot(line_cX,line_cY,color='green')
                plt.xlabel(fullmat.columns[0])
                plt.ylabel(fullmat.columns[1])
                plt.title((ccog+' '+cctrl+' '+cstate+' | SP r = '+str(corr_r)),loc='right')

                #Save.
                corrtype = 'spearman'
                outfile = (outpath+'/quad_'+ccog+'_'+cctrl+'_'+cstate+'_'+
                        cstatetype+'_'+septypeint+'_'+corrtype+'_covha_'+outlab+'.jpg')
                plt.savefig(outfile,bbox_inches='tight',dpi=360)  
                plt.close()

            #Add.
            pr_corrsave_all.append(pr_corrsave)
            sp_corrsave_all.append(sp_corrsave)
    
    #Format and save.
    corrtype = 'pearson'
    corrsave_out = pd.DataFrame(pr_corrsave_all,index=corrsave_labs,columns=csplits)
    outfile = (outpath+'/'+cstatetype+'_'+septypeint+'_'+corrtype+'_covha_corr.csv')
    corrsave_out.to_csv(outfile,index=True,header=True)
    corrtype = 'spearman'
    corrsave_out = pd.DataFrame(sp_corrsave_all,index=corrsave_labs,columns=csplits)
    outfile = (outpath+'/'+cstatetype+'_'+septypeint+'_'+corrtype+'_covha_corr.csv')
    corrsave_out.to_csv(outfile,index=True,header=True)
    print('Exploration done.')
