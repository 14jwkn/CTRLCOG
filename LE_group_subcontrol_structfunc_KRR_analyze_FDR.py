# -*- coding: utf-8 -*-
"""
For a given k for clustering, type of SC normalization, type and percentage of 
thresholding, number of CV repetitions, CV inner k, CV outer k, and number of permutations,
FDR correct across model R2 and model R2 comparisons, FDR correct across model
Haufe scores, and FDR correct across brain map Spearman's correlations, quadratic
regression R2, and Dice coefficients.

Output:
'testacc_FDR.csv' FDR-corrected p-values for model R2.
ccomp+'_FDR.csv' FDR-corrected p-values for model R2 comparison.
ccog+'_covha_featscores_FDR.csv' FDR-corrected p-values for model Haufe scores.
wantpred+'_'+cctrl+'_'+cstatetype+'_'+cseptype+'_covha_'+cblock+'_'+ccog+'_summary.csv' FDR-corrected p-values for Spearman's correlations and R2 for quadratic regression between brain maps and Haufe scores. 
'/dc_'+cctrl+'_'+cstatetype+'_'+ccog+'_covha_summary.csv' FDR-corrected p-values for Dice coefficient between brain maps and Haufe scores. 

Usage: 
    LE_group_subcontrol_structfunc_KRR_analyze_FDR.py <k> <sctype> <threstype> <thresval> <nrep> <inner_k> <outer_k> <nperm>
Arguments:
    
    <k> k for k-means clustering
    <sctype> SC type of normalization
    <threstype> Threshold type
    <thresval> Threshold value
    <nrep> Number of CV repetitions
    <inner_k> Inner K in K-fold CV hyperparameter search
    <outer_k> Outer K in K-fold CV
    <nperm> Number of permutations

"""
import os
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
from docopt import docopt

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
    nperm = args['<nperm>']
    print('Doing:',k,sctype,threstype,
          thresval,nrep,inner_k,outer_k,nperm)

    # Set base path and parameters.
    subgroup = 'full'
    sc_subgroup = 'dr_full'
    basepath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
        subgroup + '/' + k +'/SC_dFC/' + sc_subgroup + '/collect/'+
        threstype + '/' +thresval + '/')
    nk = int(k)
    
    # Read overall prediction p-values.
    ctrl_list = ['ave','mod','abs_deg']
    nctrl = len(ctrl_list)
    state_list = ['SC','dFCcat','SC_dFCcat']
    nstate = len(state_list)
    septype = 'comCFAng'
    coglist = ['gCFA','P24_CR','PV']
    ncog = len(coglist)
    raw_p = pd.Series([])
    for cctrl in ctrl_list:
        for cstate in state_list:
            inpath = (basepath +
                'KRRXFS/' +
                cctrl + '_' + cstate + '_' + septype + '_' +
                nrep + '_' + inner_k + '_' + outer_k + '_' + sctype + '/')
            infile = (inpath+'testacc.csv')
            inp = pd.read_csv(infile,index_col=0).loc[coglist,'OneP']
            inp.index = [(cctrl+' '+cstate+' '+x) for x in coglist]
            raw_p = pd.concat([raw_p,inp])

    # Read overall prediction comparison p-values.
    comp_list = ['dFC_vs_SC','SC_dFC_vs_dFC','sFC_vs_SC','dFC_vs_sFC','SC_dFC_vs_SC_sFC','cogcompare','strcompare']
    for ccomp in comp_list:
        cctrl = 'controlcomp'
        cstate = 'statecomp'
        septype = 'cogcomp'
        inpath = (basepath +
                'KRRXF/' +
                cctrl + '_' + cstate + '_' + septype + '_' +
                nrep + '_' + inner_k + '_' + outer_k + '_' + sctype + '/')
        infile = (inpath+ccomp+'.csv')
        inp = pd.read_csv(infile,index_col=0).loc[:,'P']
        inp.index = [(ccomp+' '+x) for x in inp.index]
        raw_p = pd.concat([raw_p,inp])
    
    # FDR correct.
    form_p = np.asarray(raw_p,dtype=float)
    fdr_p = false_discovery_control(form_p,method='bh')
    fdr_p = pd.Series(fdr_p,index=raw_p.index)

    # Output overall prediction p-values.
    ctrl_list = ['ave','mod','abs_deg']
    nctrl = len(ctrl_list)
    state_list = ['SC','dFCcat','SC_dFCcat']
    nstate = len(state_list)
    septype = 'comCFAng'
    coglist = ['gCFA','P24_CR','PV']
    ncog = len(coglist)
    for cctrl in ctrl_list:
        for cstate in state_list:
            inpath = (basepath +
                'KRRXFS/' +
                cctrl + '_' + cstate + '_' + septype + '_' +
                nrep + '_' + inner_k + '_' + outer_k + '_' + sctype + '/')
            infile = (inpath+'testacc.csv')
            inp = pd.read_csv(infile,index_col=0).loc[coglist,:]
            inp.index = [(cctrl+' '+cstate+' '+x) for x in coglist]
            inp.loc[:,'OneP_FDR'] = fdr_p[inp.index]
            outfile = (inpath+'testacc_FDR.csv')
            inp.to_csv(outfile)

    # Output overall prediction comparison p-values.
    comp_list = ['dFC_vs_SC','SC_dFC_vs_dFC','sFC_vs_SC','dFC_vs_sFC','SC_dFC_vs_SC_sFC','cogcompare','strcompare']
    for ccomp in comp_list:
        cctrl = 'controlcomp'
        cstate = 'statecomp'
        septype = 'cogcomp'
        inpath = (basepath +
                'KRRXF/' +
                cctrl + '_' + cstate + '_' + septype + '_' +
                nrep + '_' + inner_k + '_' + outer_k + '_' + sctype + '/')
        infile = (inpath+ccomp+'.csv')
        inp = pd.read_csv(infile,index_col=0)
        inp.index = [(ccomp+' '+x) for x in inp.index]
        inp.loc[:,'P_FDR'] = fdr_p[inp.index]
        outfile = (inpath+ccomp+'_FDR.csv')
        inp.to_csv(outfile)
    
    # Read feature importance p-values.
    ctrl_list = ['ave','mod','abs_deg']
    nctrl = len(ctrl_list)
    septype = 'comCFAng'
    cstate = 'SC_dFCcat'
    cog_list = ['gCFA','P24_CR','PV']
    ncog = len(cog_list)
    raw_p = pd.Series([])
    for cctrl in ctrl_list:
        for ccog in cog_list:
            inpath = (basepath+'KRRXFS/' +cctrl + '_' + cstate + '_' + septype + '_' +
                nrep + '_' + inner_k + '_' + outer_k + '_' + sctype + '/')
            infile = (inpath+ccog+'_covha_featscores.csv')
            inp = pd.read_csv(infile,index_col=0).loc[:,'TwoP']
            inp.index = [(cctrl+' '+ccog+' '+x) for x in inp.index]
            raw_p = pd.concat([raw_p,inp])
    
    # FDR correct.
    form_p = np.asarray(raw_p,dtype=float)
    fdr_p = false_discovery_control(form_p,method='bh')
    fdr_p = pd.Series(fdr_p,index=raw_p.index)

    # Output feature importance p-values.
    ctrl_list = ['ave','mod','abs_deg']
    nctrl = len(ctrl_list)
    septype = 'comCFAng'
    cstate = 'SC_dFCcat'
    cog_list = ['gCFA','P24_CR','PV']
    ncog = len(cog_list)
    for cctrl in ctrl_list:
        for ccog in cog_list:
            inpath = (basepath+'KRRXFS/' +cctrl + '_' + cstate + '_' + septype + '_' +
                nrep + '_' + inner_k + '_' + outer_k + '_' + sctype + '/')
            infile = (inpath+ccog+'_covha_featscores.csv')
            inp = pd.read_csv(infile,index_col=0).loc[:,['Score','TwoP']]
            raw_index = inp.index.copy()
            inp.index = [(cctrl+' '+ccog+' '+x) for x in raw_index]
            inp.loc[:,'TwoP_FDR'] = fdr_p[inp.index]
            inp.index = raw_index
            outfile = (inpath+ccog+'_covha_featscores_FDR.csv')
            inp.to_csv(outfile)

    # Read feature importance spin test p-values for quadratic regression R2 and Spearman's correlation.
    block_list = ['ctrlpred','grad1pred','areaDC']
    area_list = ['PFIT','Extended']
    xscale = 'no'
    xcon = 'no'
    pselect = 'TwoP_FDR'
    wantpred = 'full'
    wantperm = '1000'
    cstatetype = 'SC_dFCcat'
    cstatetype_list = ['sc'] + ['s'+str(x+1) for x in range(nk)]
    cseptype = 'comCFAng'
    explorepath =  (basepath+'KRRXFS/controlcomp_statecomp_cogcomp_'+
                     nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/correxplore/'+
                     xscale+'_'+xcon+'_'+pselect+'/specialp/')
    raw_p = pd.Series([])
    for cblock in block_list:
        inpath = (explorepath+cblock+'/')
        if cblock == 'ctrlpred':
            coglist = ['PV','P24_CR','gCFA']
            for ccog in coglist:
                if ccog == 'PV':
                    ctrl_list = ['ave']
                elif ccog == 'P24_CR':
                    ctrl_list = ['mod']
                elif ccog == 'gCFA':
                    ctrl_list = ['ave','mod','abs_deg']
                for cctrl in ctrl_list:
                    corrtype = 'quadreg'
                    cstat = 'r2'
                    infile = (inpath+wantpred+'/'+wantpred+'_'+cctrl+'_'+cstatetype+
                              '_'+cseptype+'_'+corrtype+'_covha_'+cblock+'_'+cstat+'_pval.csv')
                    inp = pd.read_csv(infile,index_col=0).loc[ccog,:]
                    inp.index = [(corrtype+' '+cblock+' '+ccog+' '+cctrl+' '+x) for x in inp.index]
                    raw_p = pd.concat([raw_p,inp])
                    corrtype = 'spearman'
                    infile = (inpath+wantpred+'/'+wantpred+'_'+cctrl+'_'+cstatetype+
                              '_'+cseptype+'_'+corrtype+'_covha_'+cblock+'_pval.csv')
                    inp = pd.read_csv(infile,index_col=0).loc[ccog,:]
                    inp.index = [(corrtype+' '+cblock+' '+ccog+' '+cctrl+' '+x) for x in inp.index]
                    raw_p = pd.concat([raw_p,inp]) 
        elif cblock == 'grad1pred':
            ccog = 'gCFA'
            ctrl_list = ['ave','mod','abs_deg']
            for cctrl in ctrl_list:
                corrtype = 'quadreg'
                cstat = 'r2'
                infile = (inpath+wantpred+'/'+wantpred+'_'+cctrl+'_'+cstatetype+
                            '_'+cseptype+'_'+corrtype+'_covha_'+cblock+'_'+cstat+'_pval.csv')
                inp = pd.read_csv(infile,index_col=0).loc[ccog,:]
                inp.index = [(corrtype+' '+cblock+' '+ccog+' '+cctrl+' '+x) for x in inp.index]
                raw_p = pd.concat([raw_p,inp])
                corrtype = 'spearman'
                infile = (inpath+wantpred+'/'+wantpred+'_'+cctrl+'_'+cstatetype+
                            '_'+cseptype+'_'+corrtype+'_covha_'+cblock+'_pval.csv')
                inp = pd.read_csv(infile,index_col=0).loc[ccog,:]
                inp.index = [(corrtype+' '+cblock+' '+ccog+' '+cctrl+' '+x) for x in inp.index]
                raw_p = pd.concat([raw_p,inp])
        elif cblock == 'areaDC':
            ccog = 'gCFA'
            ctrl_list = ['ave','mod','abs_deg']
            for cctrl in ctrl_list:
                infile = (inpath+'/dc_'+cctrl+'_'+cstatetype+
                            '_'+ccog+'_covha_pval.csv')
                inp_both = pd.read_csv(infile,index_col=0).loc[cstatetype_list,area_list]
                for carea in area_list:
                    inp = inp_both.loc[:,carea]
                    inp.index = [(carea+' '+cblock+' '+ccog+' '+cctrl+' '+x) for x in inp.index]
                    raw_p = pd.concat([raw_p,inp])

    # FDR correct.
    form_p = np.asarray(raw_p,dtype=float)
    fdr_p = false_discovery_control(form_p,method='bh')
    fdr_p = pd.Series(fdr_p,index=raw_p.index)

    # Output feature importance correspondence p-values.
    for cblock in block_list:
        inpath = (explorepath+cblock+'/')
        if cblock == 'ctrlpred':
            coglist = ['PV','P24_CR','gCFA']
            for ccog in coglist:
                if ccog == 'PV':
                    ctrl_list = ['ave']
                elif ccog == 'P24_CR':
                    ctrl_list = ['mod']
                elif ccog == 'gCFA':
                    ctrl_list = ['ave','mod','abs_deg']
                for cctrl in ctrl_list:

                    # Read quadratic regression R2, B2, and R2 FDR p.
                    corrtype = 'quadreg'
                    infile = (inpath+wantpred+'/'+wantpred+'_'+cctrl+'_'+cstatetype+
                              '_'+cseptype+'_'+corrtype+'_covha_'+cblock+'_r2.csv')
                    cr2 = pd.read_csv(infile,index_col=0).loc[ccog,:]
                    infile = (inpath+wantpred+'/'+wantpred+'_'+cctrl+'_'+cstatetype+
                              '_'+cseptype+'_'+corrtype+'_covha_'+cblock+'_b2.csv')
                    cb2 = pd.read_csv(infile,index_col=0).loc[ccog,:]
                    old_index = cr2.index.copy()
                    want_index = [(corrtype+' '+cblock+' '+ccog+' '+cctrl+' '+x) for x in old_index]
                    cr2_p = fdr_p[want_index]

                    # Read Spearman's correlation and FDR p.
                    corrtype = 'spearman'
                    infile = (inpath+wantpred+'/'+wantpred+'_'+cctrl+'_'+cstatetype+
                              '_'+cseptype+'_'+corrtype+'_covha_'+cblock+'.csv')
                    csp =  pd.read_csv(infile,index_col=0).loc[ccog,:]
                    want_index = [(corrtype+' '+cblock+' '+ccog+' '+cctrl+' '+x) for x in old_index]
                    csp_p = fdr_p[want_index]

                    # Reset indices and concatenate them all to save.
                    cr2.index = old_index
                    cr2_p.index = old_index
                    cb2.index = old_index
                    csp.index = old_index
                    csp_p.index = old_index
                    outcollect = pd.concat([cr2,cr2_p,cb2,csp,csp_p],axis=1)
                    outcollect.columns = ['R2','R2_P','B2','SP','SP_P']
                    outfile = (inpath+wantpred+'/'+wantpred+'_'+cctrl+'_'+cstatetype+
                              '_'+cseptype+'_covha_'+cblock+'_'+ccog+'_summary.csv')
                    outcollect.to_csv(outfile)
        elif cblock == 'grad1pred':
            ccog = 'gCFA'
            ctrl_list = ['ave','mod','abs_deg']
            for cctrl in ctrl_list:
                
                # Read quadratic regression R2, B2, and R2 FDR p.
                corrtype = 'quadreg'
                infile = (inpath+wantpred+'/'+wantpred+'_'+cctrl+'_'+cstatetype+
                            '_'+cseptype+'_'+corrtype+'_covha_'+cblock+'_r2.csv')
                cr2 = pd.read_csv(infile,index_col=0).loc[ccog,:]
                infile = (inpath+wantpred+'/'+wantpred+'_'+cctrl+'_'+cstatetype+
                            '_'+cseptype+'_'+corrtype+'_covha_'+cblock+'_b2.csv')
                cb2 = pd.read_csv(infile,index_col=0).loc[ccog,:]
                old_index = cr2.index.copy()
                want_index = [(corrtype+' '+cblock+' '+ccog+' '+cctrl+' '+x) for x in old_index]
                cr2_p = fdr_p[want_index]

                # Read Spearman's correlation and FDR p.
                corrtype = 'spearman'
                infile = (inpath+wantpred+'/'+wantpred+'_'+cctrl+'_'+cstatetype+
                            '_'+cseptype+'_'+corrtype+'_covha_'+cblock+'.csv')
                csp =  pd.read_csv(infile,index_col=0).loc[ccog,:]
                want_index = [(corrtype+' '+cblock+' '+ccog+' '+cctrl+' '+x) for x in old_index]
                csp_p = fdr_p[want_index]

                # Reset indices and concatenate them all to save.
                cr2.index = old_index
                cr2_p.index = old_index
                cb2.index = old_index
                csp.index = old_index
                csp_p.index = old_index
                outcollect = pd.concat([cr2,cr2_p,cb2,csp,csp_p],axis=1)
                outcollect.columns = ['R2','R2_P','B2','SP','SP_P']
                outfile = (inpath+wantpred+'/'+wantpred+'_'+cctrl+'_'+cstatetype+
                            '_'+cseptype+'_covha_'+cblock+'_'+ccog+'_summary.csv')
                outcollect.to_csv(outfile)
        elif cblock == 'areaDC':
            ccog = 'gCFA'
            ctrl_list = ['ave','mod','abs_deg']
            for cctrl in ctrl_list:

                # Read DC for each area, append FDR p, and save.
                infile = (inpath+'/dc_'+cctrl+'_'+cstatetype+
                            '_'+ccog+'_covha.csv')
                outcollect = pd.read_csv(infile,index_col=0).loc[cstatetype_list,area_list]
                old_index = outcollect.index.copy()
                for carea in area_list:
                    want_index = [(carea+' '+cblock+' '+ccog+' '+cctrl+' '+x) for x in old_index]
                    carea_p = fdr_p[want_index]
                    carea_p.index = old_index
                    outcols = outcollect.columns.tolist() + [(carea+'_P')]
                    outcollect = pd.concat([outcollect,carea_p],axis=1)
                    outcollect.columns = outcols
                outfile = (inpath+'/dc_'+cctrl+'_'+cstatetype+
                            '_'+ccog+'_covha_summary.csv')
                outcollect.to_csv(outfile)
