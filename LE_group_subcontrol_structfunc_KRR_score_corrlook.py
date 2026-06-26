# -*- coding: utf-8 -*-
"""
For a given k, SC type of normalization, threshold type, threshold value, and
number of repetitions, inner loop number, and outer loop number for CV, investigate
correlations to confirm correspondence - these include the correlation between
group-averaged controllability and strength, individual SC and dFC regional importance
to prediction and the combined SC & dFC regional importance for the main model, 
and controllability and strength regional importance for the main model. 
Output:
group_ctrl_corr.csv Correlation between group-averaged controllability and strength for SC, sFC, and dFC.
SC_vs_dFC_feature.csv Correlation between individual SC and dFC models regional importance and combined SC & dFC model regional importance for the main model.
ctrlcomp_feature.csv Correlation between controllability and strength regional importance for the main model.

Usage: 
    LE_group_subcontrol_structfunc_KRR_score_corrlook.py <k> <sctype> <threstype> <thresval> <nrep> <inner_k> <outer_k> 
    
Arguments:
    
    <k> k for k-means clustering
    <sctype> SC type of normalization
    <threstype> Threshold type
    <thresval> Threshold value
    <nrep> Number of CV repetitions
    <inner_k> Inner K in K-fold CV hyperparameter search
    <outer_k> Outer K in K-fold CV

"""

import os, sys, time, random
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from docopt import docopt

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
                                '_'+nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+
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
    print('Doing:',k,sctype,threstype,thresval,nrep,inner_k,outer_k)

    #Set paths.
    subgroup = 'full'
    sc_subgroup = 'dr_full'
    basepath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                subgroup+'/'+k+'/SC_dFC/'+
                sc_subgroup+'/collect/'+threstype+'/'+thresval+'/')
    scpath = ('../outputs/d_SC/'+sc_subgroup+'/'+threstype+'/'+thresval+'/'+sctype+'/')
    sFCpath = ('../outputs/r_sFC/'+sc_subgroup+'/'+threstype+'/'+thresval+'/')
    explorepath =  (basepath+'KRRXFS/controlcomp_statecomp_cogcomp_'+
                     nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/corrlook/')
    os.makedirs(explorepath,exist_ok=True)
    
    #Initialize basic values.
    nroi = 360
    nk = int(k)
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

    # Average across people for each state, then correlate between AC and MC, AC and S, and MC and S.
    ctrl_list = ['ave','mod','abs_deg']
    cstatetype = 'SC_sFC_dFCcat'
    csplits = []
    if 'SC' in cstatetype:
        csplits = csplits + ['sc']
    if 'sFC' in cstatetype:
        csplits = csplits + ['sFC']
    if 'dFC' in cstatetype:
        csplits = csplits + [('s'+str(x+1)) for x in range(nk)]
    nsplit = len(csplits)
    comp_list = ['ave-mod','ave-abs_deg','mod-abs_deg']
    ncomp = len(comp_list)
    comp_mat = pd.DataFrame(np.zeros((nsplit,ncomp)),index=csplits,columns=comp_list)
    for ccomp in comp_list:
        cctrl1, cctrl2 = ccomp.split('-')
        for cstate in csplits:
            inkey = (cctrl1 + '.' + cstate)
            cmat1 = ctrl_collect[inkey].mean()
            inkey = (cctrl2 + '.' + cstate)
            cmat2 = ctrl_collect[inkey].mean()
            comp_mat.loc[cstate,ccomp] = pearsonr(cmat1,cmat2).statistic
    
    # Save.
    outfile = (explorepath+'group_ctrl_corr.csv')
    comp_mat.to_csv(outfile)

    # Correlate between SC & dFC and SC and dFC feature scores.
    csplits = [('s'+str(x+1)) for x in range(nk)]
    nsplit = len(csplits)
    ccog = 'gCFA'
    comp_list = ['SC','dFC']
    ncomp = len(comp_list)
    ctrl_list = ['ave','mod','abs_deg']
    nctrl = len(ctrl_list)
    comp_mat = pd.DataFrame(np.zeros((nctrl,ncomp)),index=ctrl_list,columns=comp_list)
    for cctrl in ctrl_list:

        # Read in SC and dFC feature scores from SC & dFC model.
        cstatetype = 'SC_dFCcat'
        cstate = 'sc'
        inkey = (cctrl + '.' + cstatetype + '.' + ccog + '.covha.' + cstate)
        sc_both = pred_collect[inkey]
        dFC_both = []
        for cstate in csplits:
            inkey = (cctrl + '.' + cstatetype + '.' + ccog + '.covha.' + cstate)
            dFC_both.append(pred_collect[inkey])
        dFC_both = pd.concat(dFC_both)

        # Read in SC feature scores from separate model.
        cstatetype = 'SC'
        cstate = 'sc'
        inkey = (cctrl + '.' + cstatetype + '.' + ccog + '.covha.' + cstate)
        sc_sep = pred_collect[inkey]

        # Read in dFC feature scores from separate model.
        cstatetype = 'dFCcat'
        dFC_sep = []
        for cstate in csplits:
            inkey = (cctrl + '.' + cstatetype + '.' + ccog + '.covha.' + cstate)
            dFC_sep.append(pred_collect[inkey])
        dFC_sep = pd.concat(dFC_sep)

        # Correlate and save.
        comp_mat.loc[cctrl,'SC'] = pearsonr(sc_both,sc_sep).statistic
        comp_mat.loc[cctrl,'dFC'] = pearsonr(dFC_both,dFC_sep).statistic
    
    # Save.
    outfile = (explorepath+'SC_vs_dFC_feature.csv')
    comp_mat.to_csv(outfile)

    # Correlate between each control version for g SC & dFC feature scores.
    cstatetype = 'SC_dFCcat'
    csplits = []
    if 'SC' in cstatetype:
        csplits = csplits + ['sc']
    if 'sFC' in cstatetype:
        csplits = csplits + ['sFC']
    if 'dFC' in cstatetype:
        csplits = csplits + [('s'+str(x+1)) for x in range(nk)]
    nsplit = len(csplits)
    cog_list = ['gCFA','P24_CR','PV']
    ncog = len(cog_list)
    comp_list = ['ave-mod','ave-abs_deg','mod-abs_deg']
    ncomp = len(comp_list)
    comp_mat = pd.DataFrame(np.zeros((ncog,ncomp)),index=cog_list,columns=comp_list)
    for ccog in cog_list:
        for ccomp in comp_list:
            cctrl1, cctrl2 = ccomp.split('-')

            # Read.
            cfeat1 = []
            cfeat2 = []
            for cstate in csplits:
                inkey = (cctrl1 + '.' + cstatetype + '.' + ccog + '.covha.' + cstate)
                cfeat1.append(pred_collect[inkey])
                inkey = (cctrl2 + '.' + cstatetype + '.' + ccog + '.covha.' + cstate)
                cfeat2.append(pred_collect[inkey])
            cfeat1 = pd.concat(cfeat1)
            cfeat2 = pd.concat(cfeat2)
            
            # Correlate and save.
            comp_mat.loc[ccog,ccomp] = pearsonr(cfeat1,cfeat2).statistic
    
    # Save.
    outfile = (explorepath+'ctrlcomp_feature.csv')
    comp_mat.to_csv(outfile)
    print('Exploration done.')
