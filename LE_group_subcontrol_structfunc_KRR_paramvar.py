# -*- coding: utf-8 -*-
"""
For a given number of repetitions, inner loop number, and outer loop number for CV, 
collect the CV average test set accuracy score from the models using different variations 
of average or modal controllability, matrix modalities, thresholding types, thresholding values, 
and cognitive measures. First build the labels, then using the labels, 
read them all into a matrix for plotting use.
Output:
outeracc_allver.csv CV average test set accuracy scores from models using different parameters

Usage: 
    LE_group_subcontrol_structfunc_KRR_paramvar.py <nrep> <inner_k> <outer_k>
    
Arguments:
    
    <nrep> Number of CV repetitions
    <inner_k> Inner K in K-fold CV hyperparameter search
    <outer_k> Outer K in K-fold CV

"""

import os, sys, time, random
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from docopt import docopt

if __name__ == '__main__':
    __spec__ = None
    
    #Catches arguments.
    args = docopt(__doc__)
    nrep = args['<nrep>']
    inner_k = args['<inner_k>']
    outer_k = args['<outer_k>']
    nperm = args['<nperm>']
    print('Doing:',nrep,inner_k,outer_k)

    #Set paths.
    subgroup = 'full'
    sc_subgroup = 'dr_full'
    outpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
               subgroup+'/multi_param/'+sc_subgroup+'/KRRXFS/'+
               nrep+'_'+inner_k+'_'+outer_k+'/')
    os.makedirs(outpath,exist_ok=True)

    #Set up parameters.
    ks = [str(x) for x in range(2,13)]
    nk = len(ks)
    ctrltypes = ['ave','mod']
    nctrl = len(ctrltypes)
    threstypes = ['groupconsist','groupstr']
    nthrestypes = len(threstypes)
    thresvals = ['10','20','30','40','50','60','70','80','90']
    nthresvals = len(thresvals)
    septype_list = ['comCFAng',
                    'comPCAng',
                    'comEFAng',
                    'comNT']
    nseptype = len(septype_list)
    cog_septype = [['gCFA','P24_CR','PV','gFngCFA','gCngCFA'],
                   ['gPCA2','gFngPCA2','gCngPCA2'],
                   ['gEFA','gFngEFA','gCngEFA'],
                   ['NT','NF','NC','NFnNT','NCnNT']]
    coglist = sum(cog_septype,[])
    ncog = len(coglist)

    #Generate row labels for fpt.
    row_labs = []
    sctype = 'fpt'
    singstates = ['SC']
    nsing = len(singstates)
    multstates = ['dFCcat','SC_dFCcat']
    nmult = len(multstates)
    statetypes = singstates + multstates
    nstatetypes = len(statetypes)
    for ctidx in range(nctrl):
        cctrl = ctrltypes[ctidx]
        for stidx in range(nstatetypes):
            cstatetype = statetypes[stidx]

            #If k is a factor.
            if ((cstatetype == 'SC_sFC_dFCcat') or
                (cstatetype == 'SC_dFCcat') or
                (cstatetype == 'dFCcat')):

                #No thresholding.
                for kidx in range(nk):
                    ck = ks[kidx]

                    #Whether there is SC or not.
                    if cstatetype == 'dFCcat':
                        clab = ('none.'+cctrl+'.'+cstatetype+'.'+ck+'.none.0')
                    else:
                        clab = (sctype+'.'+cctrl+'.'+cstatetype+'.'+ck+'.none.0')
                    row_labs.append(clab)

                    #Each thresholding.
                    for tyidx in range(nthrestypes):
                        cthrestype = threstypes[tyidx]
                        for tvidx in range(nthresvals):
                            cthresval = thresvals[tvidx]

                            #Whether there is SC or not.
                            if cstatetype == 'dFCcat':
                                clab = ('none.'+cctrl+'.'+cstatetype+'.'+ck+'.'+cthrestype+'.'+cthresval)
                            else:
                                clab = (sctype+'.'+cctrl+'.'+cstatetype+'.'+ck+'.'+cthrestype+'.'+cthresval)
                            row_labs.append(clab)
            
            #If k is not a factor.
            else:

                #No thresholding. Whether there is SC or not.
                if cstatetype == 'dFCcat':
                    clab = ('none.'+cctrl+'.'+cstatetype+'.0.none.0')
                else:
                    clab = (sctype+'.'+cctrl+'.'+cstatetype+'.0.none.0')
                row_labs.append(clab)

                #Each thresholding.
                for tyidx in range(nthrestypes):
                    cthrestype = threstypes[tyidx]
                    for tvidx in range(nthresvals):
                        cthresval = thresvals[tvidx]

                        #Whether there is SC or not.
                        if cstatetype == 'dFCcat':
                            clab = ('none.'+cctrl+'.'+cstatetype+'.0.'+cthrestype+'.'+cthresval)
                        else:
                            clab = (sctype+'.'+cctrl+'.'+cstatetype+'.0.'+cthrestype+'.'+cthresval)
                        row_labs.append(clab)

    #Generate row labels for raw.
    sctype = 'raw'
    singstates = ['SC']
    nsing = len(singstates)
    multstates = ['SC_dFCcat']
    nmult = len(multstates)
    statetypes = singstates + multstates
    nstatetypes = len(statetypes)
    for ctidx in range(nctrl):
        cctrl = ctrltypes[ctidx]
        for stidx in range(nstatetypes):
            cstatetype = statetypes[stidx]

            #If k is a factor.
            if ((cstatetype == 'SC_sFC_dFCcat') or
                (cstatetype == 'SC_dFCcat') or
                (cstatetype == 'dFCcat')):

                #No thresholding.
                for kidx in range(nk):
                    ck = ks[kidx]

                    #Whether there is SC or not.
                    if cstatetype == 'dFCcat':
                        clab = ('none.'+cctrl+'.'+cstatetype+'.'+ck+'.none.0')
                    else:
                        clab = (sctype+'.'+cctrl+'.'+cstatetype+'.'+ck+'.none.0')
                    row_labs.append(clab)

                    #Each thresholding.
                    for tyidx in range(nthrestypes):
                        cthrestype = threstypes[tyidx]
                        for tvidx in range(nthresvals):
                            cthresval = thresvals[tvidx]

                            #Whether there is SC or not.
                            if cstatetype == 'dFCcat':
                                clab = ('none.'+cctrl+'.'+cstatetype+'.'+ck+'.'+cthrestype+'.'+cthresval)
                            else:
                                clab = (sctype+'.'+cctrl+'.'+cstatetype+'.'+ck+'.'+cthrestype+'.'+cthresval)
                            row_labs.append(clab)
            
            #If k is not a factor.
            else:

                #No thresholding. Whether there is SC or not.
                if cstatetype == 'dFCcat':
                    clab = ('none.'+cctrl+'.'+cstatetype+'.0.none.0')
                else:
                    clab = (sctype+'.'+cctrl+'.'+cstatetype+'.0.none.0')
                row_labs.append(clab)

                #Each thresholding.
                for tyidx in range(nthrestypes):
                    cthrestype = threstypes[tyidx]
                    for tvidx in range(nthresvals):
                        cthresval = thresvals[tvidx]

                        #Whether there is SC or not.
                        if cstatetype == 'dFCcat':
                            clab = ('none.'+cctrl+'.'+cstatetype+'.0.'+cthrestype+'.'+cthresval)
                        else:
                            clab = (sctype+'.'+cctrl+'.'+cstatetype+'.0.'+cthrestype+'.'+cthresval)
                        row_labs.append(clab)

    #Gather the average outer accuracy values of each hyperparameter set for 
    #all cognitive variables.
    inkey = ('/rep_testacc')
    print('Averaged outer accuracy.')
    cacc = pd.DataFrame(np.zeros((len(row_labs),ncog)),index=row_labs,columns=coglist)

    #Do fpt with the same labels.
    sctype = 'fpt'
    singstates = ['SC']
    nsing = len(singstates)
    multstates = ['dFCcat','SC_dFCcat']
    nmult = len(multstates)
    statetypes = singstates + multstates
    nstatetypes = len(statetypes)
    for ctidx in range(nctrl):
        cctrl = ctrltypes[ctidx]
        for stidx in range(nstatetypes):
            cstatetype = statetypes[stidx]

            #If k is a factor.
            if ((cstatetype == 'SC_sFC_dFCcat') or
                (cstatetype == 'SC_dFCcat') or
                (cstatetype == 'dFCcat')):

                #No thresholding.
                for kidx in range(nk):
                    ck = ks[kidx]

                    #Whether there is SC or not.
                    if cstatetype == 'dFCcat':
                        clab = ('none.'+cctrl+'.'+cstatetype+'.'+ck+'.none.0')
                    else:
                        clab = (sctype+'.'+cctrl+'.'+cstatetype+'.'+ck+'.none.0')
                    
                    #For each cognitive set, read and place in.
                    for cseptype in septype_list:
                        inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                                subgroup+'/'+ck+'/SC_dFC/'+
                                sc_subgroup+'/collect/none/0/KRRXFS/'+
                                cctrl+'_'+cstatetype+'_'+cseptype+
                                '_'+nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/')
                        infile = (inpath+'score_collect.h5')
                        store = pd.HDFStore(infile,'r')
                        testacc = store.select(inkey)
                        store.close()
                        colnames = testacc.columns
                        cacc.loc[clab,colnames] = testacc.mean(axis=0)
                    
                    #Each thresholding.
                    for tyidx in range(nthrestypes):
                        cthrestype = threstypes[tyidx]
                        for tvidx in range(nthresvals):
                            cthresval = thresvals[tvidx]

                            #Whether there is SC or not.
                            if cstatetype == 'dFCcat':
                                clab = ('none.'+cctrl+'.'+cstatetype+'.'+ck+'.'+cthrestype+'.'+cthresval)
                            else:
                                clab = (sctype+'.'+cctrl+'.'+cstatetype+'.'+ck+'.'+cthrestype+'.'+cthresval)

                            #For each cognitive set, read and place in.
                            for cseptype in septype_list:
                                inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                                        subgroup+'/'+ck+'/SC_dFC/'+
                                        sc_subgroup+'/collect/'+cthrestype+'/'+cthresval+'/KRRXF/'+
                                        cctrl+'_'+cstatetype+'_'+cseptype+
                                        '_'+nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/')
                                infile = (inpath+'score_collect.h5')
                                store = pd.HDFStore(infile,'r')
                                testacc = store.select(inkey)
                                store.close()
                                colnames = testacc.columns
                                cacc.loc[clab,colnames] = testacc.mean(axis=0)
 
            #If k is not a factor.
            else:

                #No thresholding.  Whether there is SC or not.
                if cstatetype == 'dFCcat':
                    clab = ('none.'+cctrl+'.'+cstatetype+'.0.none.0')
                else:
                    clab = (sctype+'.'+cctrl+'.'+cstatetype+'.0.none.0')

                #For each cognitive set, read and place in.
                for cseptype in septype_list:
                    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                            subgroup+'/6/SC_dFC/'+
                            sc_subgroup+'/collect/none/0/KRRXF/'+
                            cctrl+'_'+cstatetype+'_'+cseptype+'_'+
                            nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/')
                    infile = (inpath+'score_collect.h5')
                    store = pd.HDFStore(infile,'r')
                    testacc = store.select(inkey)
                    store.close()
                    colnames = testacc.columns
                    cacc.loc[clab,colnames] = testacc.mean(axis=0)

                #Each thresholding.
                for tyidx in range(nthrestypes):
                    cthrestype = threstypes[tyidx]
                    for tvidx in range(nthresvals):
                        cthresval = thresvals[tvidx]

                        #Whether there is SC or not.
                        if cstatetype == 'dFCcat':
                            clab = ('none.'+cctrl+'.'+cstatetype+'.0.'+cthrestype+'.'+cthresval)
                        else:
                            clab = (sctype+'.'+cctrl+'.'+cstatetype+'.0.'+cthrestype+'.'+cthresval)

                        #For each cognitive set, read and place in.
                        for cseptype in septype_list:
                            inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                                    subgroup+'/6/SC_dFC/'+
                                    sc_subgroup+'/collect/'+cthrestype+'/'+cthresval+'/KRRXF/'+
                                    cctrl+'_'+cstatetype+'_'+cseptype+'_'+
                                    nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/')
                            infile = (inpath+'score_collect.h5')
                            store = pd.HDFStore(infile,'r')
                            testacc = store.select(inkey)
                            store.close()
                            colnames = testacc.columns
                            cacc.loc[clab,colnames] = testacc.mean(axis=0)

    #Do raw with the same labels.
    sctype = 'raw'
    singstates = ['SC']
    nsing = len(singstates)
    multstates = ['SC_dFCcat']
    nmult = len(multstates)
    statetypes = singstates + multstates
    nstatetypes = len(statetypes)
    for ctidx in range(nctrl):
        cctrl = ctrltypes[ctidx]
        for stidx in range(nstatetypes):
            cstatetype = statetypes[stidx]

            #If k is a factor.
            if ((cstatetype == 'SC_sFC_dFCcat') or
                (cstatetype == 'SC_dFCcat') or
                (cstatetype == 'dFCcat')):

                #No thresholding.
                for kidx in range(nk):
                    ck = ks[kidx]

                    #Whether there is SC or not.
                    if cstatetype == 'dFCcat':
                        clab = ('none.'+cctrl+'.'+cstatetype+'.'+ck+'.none.0')
                    else:
                        clab = (sctype+'.'+cctrl+'.'+cstatetype+'.'+ck+'.none.0')
                    
                    #For each cognitive set, read and place in.
                    for cseptype in septype_list:
                        inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                                subgroup+'/'+ck+'/SC_dFC/'+
                                sc_subgroup+'/collect/none/0/KRRXF/'+
                                cctrl+'_'+cstatetype+'_'+cseptype+'_'+
                                nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/')
                        infile = (inpath+'score_collect.h5')
                        store = pd.HDFStore(infile,'r')
                        testacc = store.select(inkey)
                        store.close()
                        colnames = testacc.columns
                        cacc.loc[clab,colnames] = testacc.mean(axis=0)
                    
                    #Each thresholding.
                    for tyidx in range(nthrestypes):
                        cthrestype = threstypes[tyidx]
                        for tvidx in range(nthresvals):
                            cthresval = thresvals[tvidx]

                            #Whether there is SC or not.
                            if cstatetype == 'dFCcat':
                                clab = ('none.'+cctrl+'.'+cstatetype+'.'+ck+'.'+cthrestype+'.'+cthresval)
                            else:
                                clab = (sctype+'.'+cctrl+'.'+cstatetype+'.'+ck+'.'+cthrestype+'.'+cthresval)

                            #For each cognitive set, read and place in.
                            for cseptype in septype_list:
                                inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                                        subgroup+'/'+ck+'/SC_dFC/'+
                                        sc_subgroup+'/collect/'+cthrestype+'/'+cthresval+'/KRRXF/'+
                                        cctrl+'_'+cstatetype+'_'+cseptype+'_'+
                                        nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/')
                                infile = (inpath+'score_collect.h5')
                                store = pd.HDFStore(infile,'r')
                                testacc = store.select(inkey)
                                store.close()
                                colnames = testacc.columns
                                cacc.loc[clab,colnames] = testacc.mean(axis=0)
 
            #If k is not a factor.
            else:

                #No thresholding.  Whether there is SC or not.
                if cstatetype == 'dFCcat':
                    clab = ('none.'+cctrl+'.'+cstatetype+'.0.none.0')
                else:
                    clab = (sctype+'.'+cctrl+'.'+cstatetype+'.0.none.0')

                #For each cognitive set, read and place in.
                for cseptype in septype_list:
                    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                            subgroup+'/6/SC_dFC/'+
                            sc_subgroup+'/collect/none/0/KRRXF/'+
                            cctrl+'_'+cstatetype+'_'+cseptype+'_'+
                            nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/')
                    infile = (inpath+'score_collect.h5')
                    store = pd.HDFStore(infile,'r')
                    testacc = store.select(inkey)
                    store.close()
                    colnames = testacc.columns
                    cacc.loc[clab,colnames] = testacc.mean(axis=0)

                #Each thresholding.
                for tyidx in range(nthrestypes):
                    cthrestype = threstypes[tyidx]
                    for tvidx in range(nthresvals):
                        cthresval = thresvals[tvidx]

                        #Whether there is SC or not.
                        if cstatetype == 'dFCcat':
                            clab = ('none.'+cctrl+'.'+cstatetype+'.0.'+cthrestype+'.'+cthresval)
                        else:
                            clab = (sctype+'.'+cctrl+'.'+cstatetype+'.0.'+cthrestype+'.'+cthresval)

                        #For each cognitive set, read and place in.
                        for cseptype in septype_list:
                            inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                                    subgroup+'/6/SC_dFC/'+
                                    sc_subgroup+'/collect/'+cthrestype+'/'+cthresval+'/KRRXF/'+
                                    cctrl+'_'+cstatetype+'_'+cseptype+'_'+
                                    nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/')
                            infile = (inpath+'score_collect.h5')
                            store = pd.HDFStore(infile,'r')
                            testacc = store.select(inkey)
                            store.close()
                            colnames = testacc.columns
                            cacc.loc[clab,colnames] = testacc.mean(axis=0)                    
    
    #Package.
    cacc.loc['ValMax',:] = cacc.max(axis=0)
    cacc.loc['IDMax',:] = cacc.idxmax(axis=0)
    outfile = (outpath+'outeracc_allver.csv')
    cacc.to_csv(outfile,index=True,header=True)
    print('Pseudotuning done.')
