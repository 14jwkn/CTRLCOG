# -*- coding: utf-8 -*-
"""
For a given k for clustering, type of SC normalization, and type and percentage of thresholding, do permutation testing on the
model predicting each cognitive variable in the given batch using either average controllability, modal controllability,
or degree calculated from the set of state matrices specified. The one-sided p-value is produced from the test set CV accuracy
score because the Coefficient of Determination (COD) is only positive, and two-sided p-values are produced for the Haufe scores
because Haufe scores can be positive or negative. Accuracy, decomposed accuracy, and Haufe scores are all packaged for retrieval.
Output:
testacc.csv Test set CV accuracy score and one-sided permutation p-value.
csplitver+'.csv' Decomposed CV accuracy scores collected.
coglist[cidx]+'_covha_featscores.csv' Haufe scores and two-sided permutation p-values for each cognitive variable.

Usage: 
    LE_group_subcontrol_structfunc_KRR_analyze.py <k> <sctype> <threstype> <thresval> <controltype> <statetype> <septype> <nrep> <inner_k> <outer_k> <nperm>
    
Arguments:
    
    <k> k for k-means clustering
    <sctype> SC type of normalization
    <threstype> Threshold type
    <thresval> Threshold value
    <controltype> Control type
    <statetype> Which state 
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
from docopt import docopt

#One-sided significance.
def onep_feat(cfeat,cperm_max,cperm_min):

    #Extract.
    [nperm,nfeat] = np.shape(cperm_max)
    featlabs = cfeat.index

    #Go through columns.
    feat_pval = pd.DataFrame(np.zeros((1,nfeat)))
    feat_pval.columns = featlabs
    for fidx in range(nfeat):

        #Extract.
        featlab = featlabs[fidx]
        onetruefeat = cfeat.loc[featlab]
        onepermfeat_max = cperm_max.loc[:,featlab]
        onepermfeat_min = cperm_min.loc[:,featlab]

        #If positive or negative.
        if onetruefeat >= 0:
            permsig = ((sum(onepermfeat_max >= onetruefeat)) + 1)/(nperm + 1)   
        else:
            permsig = ((sum(onepermfeat_min <= onetruefeat)) + 1)/(nperm + 1)
            
        #Append.
        feat_pval.loc[0,featlab] = permsig
    
    #Return p-values.
    return feat_pval

#Two-sided significance.
def twop_feat(cfeat,cperm):

    #Extract.
    [nperm,nfeat] = np.shape(cperm)
    featlabs = cperm.columns

    #Absolute value each and find surpassing values.
    cfeat_abs = cfeat.abs()
    cperm_abs = cperm.abs()
    feat_pval_two = cperm_abs.ge(cfeat_abs,axis=1).sum(axis=0).add(1).div((nperm+1))
    feat_pval_two = pd.DataFrame(feat_pval_two).T
    feat_pval_two.columns = featlabs

    #Return p-values.
    return feat_pval_two

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
    krrpath = (basepath+'KRRXFS/'+controltype+'_'+statetype+'_'+septype+'_'+
               nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/')
    permpath = (krrpath+'perm_'+nperm+'/permsep/')
    
    #Initialize values.
    nroi = 360
    nk = int(k)
    inner_k = int(inner_k)
    outer_k = int(outer_k) 
    nrep = int(nrep)
    nperm = int(nperm)
    if septype == 'comCFAng':
        coglist = ['gCFA','P24_CR','PV','gFngCFA','gCngCFA']
    elif septype == 'comPCAng':
        coglist = ['gPCA2','gFngPCA2','gCngPCA2']
    elif septype == 'comEFAng':
        coglist = ['gEFA','gFngEFA','gCngEFA']
    elif septype == 'comNT':
        coglist = ['NT','NF','NC','NFnNT','NCnNT']
    ncog = len(coglist)

    #Get the space names.
    if statetype == 'SC_sFC_dFCcat':
        feature_names = ['sc','sFC','dFC']
    if statetype == 'SC_dFCcat':
        feature_names = ['sc','dFC']
    if statetype == 'SC_sFC':
        feature_names = ['sc','sFC']
    if statetype == 'dFCcat':
        feature_names = ['dFC']
    if statetype == 'SC':
        feature_names = ['sc']
    if statetype == 'sFC':
        feature_names = ['sFC']
    nspace = len(feature_names)

    #Read in feature labels.
    infile = (krrpath+'score_collect.h5')
    store = pd.HDFStore(infile,'r')
    inkey = ('/rep_covha_feat_'+coglist[0])
    featlabs = store.select(inkey).columns
    store.close()
    nfeat = len(featlabs)

    #Read in the test accuracy.
    infile = (krrpath+'score_collect.h5')
    store = pd.HDFStore(infile,'r')
    inkey = ('/rep_testacc')
    testacc = store.select(inkey)
    store.close()
    testacc = testacc.mean(axis=0)

    #Read in test accuracy permutations.
    perm_testacc = np.zeros((nperm,ncog))
    for pidx in range(int(nperm)):
        permidx = str(pidx+1)
        infile = (permpath+'cperm_'+permidx+'.h5')
        inkey = ('/rep_testacc')
        store = pd.HDFStore(infile,'r')
        accmat = store.select(inkey)
        store.close()
        accmat = accmat.mean(axis=0)
        perm_testacc[pidx,:] = accmat
    perm_testacc = pd.DataFrame(perm_testacc,columns=coglist)
    
    #Test for test accuracy significance. One-sided because COD is only positive.
    testacc_OneP = perm_testacc.ge(testacc,axis=1).sum(axis=0).add(1).div((nperm+1))

    #Package and save.
    testacc_out = pd.concat([testacc,testacc_OneP],axis=1)
    testacc_out.columns = ['Accuracy','OneP']
    outfile = (krrpath+'testacc.csv')
    testacc_out.to_csv(outfile,index=True,header=True)

    #Read in split accuracy versions.
    splitvers = ['splitacc','splitacc_strict']
    nsplitvers = len(splitvers)
    for spidx in range(nsplitvers):
        csplitver = splitvers[spidx]

        #Read in the matrix.
        infile = (krrpath+'score_collect.h5')
        store = pd.HDFStore(infile,'r')
        inkey = ('/rep_'+csplitver)
        csplitmat = store.select(inkey)
        store.close()
        csplitmat = csplitmat.mean(axis=0)
        
        #Reshape.
        csplit_out = pd.DataFrame(csplitmat.values.reshape(ncog,nspace).T,index=feature_names,columns=coglist)
        csplit_out.loc['Total'] = csplit_out.sum()

        #Save.
        outfile = (krrpath+csplitver+'.csv')
        csplit_out.to_csv(outfile,index=True,header=True)

    #Do permutation for features for each target. Do two-sided because Haufe scores can be positive and negative.
    featscores = np.zeros((nfeat,ncog))
    feat_TwoP = np.zeros((nfeat,ncog))
    for cidx in range(ncog):
        ccog = coglist[cidx]

        #Read in the feature values for this target.
        infile = (krrpath+'score_collect.h5')
        store = pd.HDFStore(infile,'r')
        inkey = ('/rep_covha_feat_'+ccog)
        cfeatvals = store.select(inkey)
        store.close()
        cfeatvals = cfeatvals.mean(axis=0)
        featscores[:,cidx] = cfeatvals

        #Read in the permutation values for this target.
        perm_featvals = np.zeros((nperm,nfeat))
        for pidx in range(int(nperm)):
            permidx = str(pidx+1)
            infile = (permpath+'cperm_'+permidx+'.h5')
            inkey = ('/rep_covha_feat_'+ccog)
            store = pd.HDFStore(infile,'r')
            cfeatperm = store.select(inkey)
            store.close()
            cfeatperm = cfeatperm.mean(axis=0)
            perm_featvals[pidx,:] = cfeatperm
        perm_featvals = pd.DataFrame(perm_featvals,columns=featlabs)

        #Find two-sided significance.
        cfeat = cfeatvals
        cperm = perm_featvals
        feat_pval = twop_feat(cfeat,cperm)
        feat_TwoP[:,cidx] = feat_pval
        print('Two-sided significant features for',ccog,': ',
            np.sum(np.transpose(feat_pval < 0.05)[0]),'of',
            np.shape(feat_pval)[1])

        #Package and save.
        for cidx in range(ncog):
            feat_out = pd.concat([pd.DataFrame(featscores[:,cidx]),
                                  pd.DataFrame(feat_TwoP[:,cidx])],axis=1)
            feat_out.index = featlabs
            feat_out.columns = ['Score','TwoP']
            outfile = (krrpath+coglist[cidx]+'_covha_featscores.csv')
            feat_out.to_csv(outfile,index=True,header=True)
    print('Analysis done.')
