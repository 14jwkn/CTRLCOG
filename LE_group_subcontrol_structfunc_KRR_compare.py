# -*- coding: utf-8 -*-
"""
For a given k for clustering, type of SC normalization, and type and percentage of thresholding, do corrected paired t-tests between
the different models predicting each cognitive variable using either average controllability, modal controllability,
or degree calculated from the set of state matrices specified. 

Output:
dFC_vs_SC.csv Compare dFC and SC models for each of the cognitive variables and controllability versions.
SC_dFC_vs_dFC Compare SC & dFC and SC models for each of the cognitive variables and controllability versions.
sFC_vs_SC.csv Compare sFC and SC models for each of the cognitive variables and controllability versions.
dFC_vs_sFC.csv Compare dFC and sFC models for each of the cognitive variables and controllability versions.
SC_dFC_vs_SC_sFC.csv Compare SC & dFC and SC & sFC models for each of the cognitive variables and controllability versions.
cogcompare.csv Compare g with gF, gC, gF-g, and gC-g for the SC & dFC model for each controllability version.
strcompare.csv Compare AC and MC with degree for the SC, dFC, and SC & dFC main models for each controllability fersion.

Usage: 
    LE_group_subcontrol_structfunc_KRR_compare.py <k> <sctype> <threstype> <thresval> <controltype> <statetype> <septype> <nrep> <inner_k> <outer_k> <nperm>
    
Arguments:
    
    <k> k for k-means clustering
    <sctype> SC type of normalization
    <threstype> Threshold type
    <thresval> Threshold value
    <controltype> Control type
    <statetype> Which states
    <septype> Cognitive variables batch
    <nrep> Number of CV repetitions
    <inner_k> Inner K in K-fold CV hyperparameter search
    <outer_k> Outer K in K-fold CV
    <nperm> Number of permutations

"""

import os, sys, time, random
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import t
from docopt import docopt

#Do corrected paired t-test for accuracy vectors. Adapted from CBIG.
def corrected_t_test(acc1,acc2):
    acc_tab = acc1 - acc2
    nullval = 0
    c_nfold = np.shape(acc_tab)[1]
    portion = (1/c_nfold)/(1-(1/c_nfold))
    c_nrep = np.shape(acc_tab)[0]
    acc_vec = acc_tab.reshape((c_nfold*c_nrep,1))
    c_nvec = np.shape(acc_vec)[0]
    corrected_var = ((1/c_nvec)+portion)*np.var(acc_vec,ddof=1)
    mu = np.mean(acc_vec)
    tval = (mu-nullval)/np.sqrt(corrected_var)
    pval = 2 * t.cdf(-np.abs(tval),(c_nvec-1))

    #Return mean of differences and p-value.
    return mu, pval

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
    nperm = args['<nperm>']
    print('Doing:',k,sctype,threstype,
          thresval,controltype,statetype,septype,nrep,
          inner_k,outer_k,nperm)

    #Set paths.
    subgroup = 'full'
    sc_subgroup = 'dr_full'
    basepath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                subgroup+'/'+k+'/SC_dFC/'+
                sc_subgroup+'/collect/'+threstype+'/'+thresval+'/')
    outpath = (basepath+'KRRXFS/'+controltype+'_'+statetype+'_'+septype+'_'+
               nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/')
    os.makedirs(outpath,exist_ok=True)
    
    #Initialize values.
    replabs = ['r'+str(i+1) for i in range(int(nrep))]
    foldlabs = ['f'+str(i+1) for i in range(int(outer_k))]

    #Set lists for the models of interest.
    statelist = ['SC','sFC','SC_sFC','dFCcat','SC_dFCcat','SC_sFC_dFCcat']
    nstate = len(statelist)
    ctrl_list = ['ave','mod','deg']
    nctrl = len(ctrl_list)
    coglist = ['gCFA','P24_CR','PV','gFngCFA','gCngCFA']
    ncog = len(coglist)
    cogset_list = ['comCFAng']
    cogset_cog = [['gCFA','P24_CR','PV','gFngCFA','gCngCFA']]
    ncogset = len(cogset_list)
 
    #Extract test accuracy matrices.
    acc_dict = {}
    for cstate in statelist:
        for cctrl in ctrl_list:
            for ccog in coglist:

                #Extract the cognitive set.
                for cogset_idx in range(ncogset):
                    ccogset_cog = cogset_cog[cogset_idx]
                    if ccog in ccogset_cog:
                        cseptype = cogset_list[cogset_idx]

                #Read in matrix.
                inpath = (basepath+'KRRXF/'+cctrl+'_'+cstate+'_'+cseptype+'_'+
                          nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/')
                infile = (inpath+'score_collect.h5')
                store = pd.HDFStore(infile,'r')
                inkey = ('/rep_testacc')
                testacc = store.select(inkey)
                store.close()
                
                #Reformat into nrep x nfold.
                ctestacc = testacc.loc[:,ccog]
                ctestacc = pd.DataFrame(ctestacc.values.reshape((int(nrep),int(outer_k))),
                                        index=replabs,columns=foldlabs)
                
                #Append.
                clab = (cstate+'_'+cctrl+'_'+ccog)
                acc_dict[clab] = ctestacc
    
    #Comparison of dFC - SC.
    cstate1 = 'dFCcat'
    cstate2 = 'SC'
    cctrl_list = ['ave','mod']
    ccog_list = ['gCFA','P24_CR','PV','gFngCFA','gCngCFA']
    row_labs = []
    for cctrl in cctrl_list:
        for ccog in ccog_list:
            row_labs.append(cctrl+' '+ccog)
    nrow = len(row_labs)
    paircom = pd.DataFrame(np.zeros((nrow,2)),index=row_labs,columns=['Mean','P'])
    for cctrl in cctrl_list:
        for ccog in ccog_list:
            crowlab = (cctrl+' '+ccog)
            
            #Extract.
            clab = (cstate1+'_'+cctrl+'_'+ccog)
            acc1 = np.array(acc_dict[clab])
            clab = (cstate2+'_'+cctrl+'_'+ccog)
            acc2 = np.array(acc_dict[clab])

            #Do test and save.
            mu,pval = corrected_t_test(acc1,acc2)
            paircom.loc[crowlab,'Mean'] = mu
            paircom.loc[crowlab,'P'] = pval
    outfile = (outpath+'dFC_vs_SC.csv')
    paircom.to_csv(outfile)

    #Comparison of SC+dFC - dFC.
    cstate1 = 'SC_dFCcat'
    cstate2 = 'dFCcat'
    cctrl_list = ['ave','mod']
    ccog_list = ['gCFA','P24_CR','PV','gFngCFA','gCngCFA']
    row_labs = []
    for cctrl in cctrl_list:
        for ccog in ccog_list:
            row_labs.append(cctrl+' '+ccog)
    nrow = len(row_labs)
    paircom = pd.DataFrame(np.zeros((nrow,2)),index=row_labs,columns=['Mean','P'])
    for cctrl in cctrl_list:
        for ccog in ccog_list:
            crowlab = (cctrl+' '+ccog)
            
            #Extract.
            clab = (cstate1+'_'+cctrl+'_'+ccog)
            acc1 = np.array(acc_dict[clab])
            clab = (cstate2+'_'+cctrl+'_'+ccog)
            acc2 = np.array(acc_dict[clab])

            #Do test and save.
            mu,pval = corrected_t_test(acc1,acc2)
            paircom.loc[crowlab,'Mean'] = mu
            paircom.loc[crowlab,'P'] = pval
    outfile = (outpath+'SC_dFC_vs_dFC.csv')
    paircom.to_csv(outfile)

    #Comparison of sFC - SC.
    cstate1 = 'sFC'
    cstate2 = 'SC'
    cctrl_list = ['ave','mod']
    ccog_list = ['gCFA','P24_CR','PV','gFngCFA','gCngCFA']
    row_labs = []
    for cctrl in cctrl_list:
        for ccog in ccog_list:
            row_labs.append(cctrl+' '+ccog)
    nrow = len(row_labs)
    paircom = pd.DataFrame(np.zeros((nrow,2)),index=row_labs,columns=['Mean','P'])
    for cctrl in cctrl_list:
        for ccog in ccog_list:
            crowlab = (cctrl+' '+ccog)
            
            #Extract.
            clab = (cstate1+'_'+cctrl+'_'+ccog)
            acc1 = np.array(acc_dict[clab])
            clab = (cstate2+'_'+cctrl+'_'+ccog)
            acc2 = np.array(acc_dict[clab])

            #Do test and save.
            mu,pval = corrected_t_test(acc1,acc2)
            paircom.loc[crowlab,'Mean'] = mu
            paircom.loc[crowlab,'P'] = pval
    outfile = (outpath+'sFC_vs_SC.csv')
    paircom.to_csv(outfile)

    #Comparison of dFC - sFC.
    cstate1 = 'dFCcat'
    cstate2 = 'sFC'
    cctrl_list = ['ave','mod']
    ccog_list = ['gCFA','P24_CR','PV','gFngCFA','gCngCFA']
    row_labs = []
    for cctrl in cctrl_list:
        for ccog in ccog_list:
            row_labs.append(cctrl+' '+ccog)
    nrow = len(row_labs)
    paircom = pd.DataFrame(np.zeros((nrow,2)),index=row_labs,columns=['Mean','P'])
    for cctrl in cctrl_list:
        for ccog in ccog_list:
            crowlab = (cctrl+' '+ccog)
            
            #Extract.
            clab = (cstate1+'_'+cctrl+'_'+ccog)
            acc1 = np.array(acc_dict[clab])
            clab = (cstate2+'_'+cctrl+'_'+ccog)
            acc2 = np.array(acc_dict[clab])

            #Do test and save.
            mu,pval = corrected_t_test(acc1,acc2)
            paircom.loc[crowlab,'Mean'] = mu
            paircom.loc[crowlab,'P'] = pval
    outfile = (outpath+'dFC_vs_sFC.csv')
    paircom.to_csv(outfile)

    #Comparison of SC+dFC - SC+sFC.
    cstate1 = 'SC_dFCcat'
    cstate2 = 'SC_sFC'
    cctrl_list = ['ave','mod']
    ccog_list = ['gCFA','P24_CR','PV','gFngCFA','gCngCFA']
    row_labs = []
    for cctrl in cctrl_list:
        for ccog in ccog_list:
            row_labs.append(cctrl+' '+ccog)
    nrow = len(row_labs)
    paircom = pd.DataFrame(np.zeros((nrow,2)),index=row_labs,columns=['Mean','P'])
    for cctrl in cctrl_list:
        for ccog in ccog_list:
            crowlab = (cctrl+' '+ccog)
            
            #Extract.
            clab = (cstate1+'_'+cctrl+'_'+ccog)
            acc1 = np.array(acc_dict[clab])
            clab = (cstate2+'_'+cctrl+'_'+ccog)
            acc2 = np.array(acc_dict[clab])

            #Do test and save.
            mu,pval = corrected_t_test(acc1,acc2)
            paircom.loc[crowlab,'Mean'] = mu
            paircom.loc[crowlab,'P'] = pval
    outfile = (outpath+'SC_dFC_vs_SC_sFC.csv')
    paircom.to_csv(outfile)

    #Comparison of g with gF, gC, gF-g and gC-g.
    cctrl_list = ['ave','mod']
    cstate_list = ['SC','dFCcat','SC_dFCcat']
    ccog_list1 = ['gCFA']
    ccog_list2 = ['P24_CR','PV','gFngCFA','gCngCFA']
    row_labs = []
    for cctrl in cctrl_list:
        for cstate in cstate_list:
            for ccog1 in ccog_list1:
                for ccog2 in ccog_list2:
                    row_labs.append(cctrl+' '+cstate+' '+ccog1+' '+ccog2)
    nrow = len(row_labs)
    paircom = pd.DataFrame(np.zeros((nrow,2)),index=row_labs,columns=['Mean','P'])
    for cctrl in cctrl_list:
        for cstate in cstate_list:
            for ccog1 in ccog_list1:
                for ccog2 in ccog_list2:
                    crowlab = (cctrl+' '+cstate+' '+ccog1+' '+ccog2)

                    #Extract.
                    clab = (cstate+'_'+cctrl+'_'+ccog1)
                    acc1 = np.array(acc_dict[clab])
                    clab = (cstate+'_'+cctrl+'_'+ccog2)
                    acc2 = np.array(acc_dict[clab])

                    #Do test and save.
                    mu,pval = corrected_t_test(acc1,acc2)
                    paircom.loc[crowlab,'Mean'] = mu
                    paircom.loc[crowlab,'P'] = pval
    outfile = (outpath+'cogcompare.csv')
    paircom.to_csv(outfile)

    #Comparison of AC and MC with degree.
    cctrl_list = ['ave','mod']
    cstate_list = ['SC','dFCcat','SC_dFCcat']
    ccoglist = ['gCFA','P24_CR','PV','gFngCFA','gCngCFA']
    row_labs = []
    for cctrl in cctrl_list:
        for cstate in cstate_list:
            for ccog in ccog_list:
                row_labs.append(cctrl+' '+cstate+' '+ccog)
    nrow = len(row_labs)
    paircom = pd.DataFrame(np.zeros((nrow,2)),index=row_labs,columns=['Mean','P'])
    for cctrl in cctrl_list:
        for cstate in cstate_list:
            for ccog in ccog_list:
                crowlab = (cctrl+' '+cstate+' '+ccog)

                #Extract.
                clab = (cstate+'_'+cctrl+'_'+ccog)
                acc1 = np.array(acc_dict[clab])
                clab = (cstate+'_deg_'+ccog)
                acc2 = np.array(acc_dict[clab])

                #Do test and save.
                mu,pval = corrected_t_test(acc1,acc2)
                paircom.loc[crowlab,'Mean'] = mu
                paircom.loc[crowlab,'P'] = pval
    
    #Save.      
    outfile = (outpath+'strcompare.csv')
    paircom.to_csv(outfile)
    print('Comparison done.')
