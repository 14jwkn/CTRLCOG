#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Relax the fMRI requirement to increase sample size to improve factor analysis. 
Get corresponding cognitive scores.
Adapts code from:
https://github.com/adolphslab/HCP_MRI-behavior/blob/master/intelligence.ipynb
Output:
cogscores.csv Contains 10 cognitive scores with expanded subjects.
gPCA_load.csv Contains the loadings from PCA.

"""

import os
import pandas as pd
import numpy as np
    
if __name__ == '__main__':
    __spec__ = None

    #Get cognitive scores.
    openpath = '../inputs/data/hcp/HCP1200_Behavioral.csv'
    df = pd.read_csv(openpath,header=0)

    #Keep subjects from all releases.
    keepSub = ((df['Release'] == 'Q1') | (df['Release'] == 'Q2') | (df['Release'] == 'Q3') 
           | (df['Release'] == 'S500') | (df['Release'] == 'S900') | (df['Release'] == 'S1200') 
           | (df['Release'] == 'MEG2'))
    print('Selected {} subjects'.format(np.sum(keepSub)))

    #Find those with no missing cognitive values.
    cogtests = ['CardSort_Unadj',
                'Flanker_Unadj',
                'ListSort_Unadj',
                'PicSeq_Unadj',
                'PicVocab_Unadj',
                'ProcSpeed_Unadj',
                'ReadEng_Unadj',
                'PMAT24_A_CR',
                'IWRD_TOT',
                'VSPLOT_TC']
    for ccog in cogtests:
        keepSub = np.logical_and(keepSub,np.logical_not(np.isnan(df[ccog])))
    print('Kept {} subjects after removing missing values'.format(np.sum(keepSub)))

    #Check mismatches with our main subject set.
    theirsub = list(map(str,df.loc[keepSub,'Subject'].values.tolist()))
    subfile = 'r_full_submain.txt'
    with open(subfile) as f:
        oursub = [subject.rstrip() for subject in f]
    mismatch = [x for x in oursub if x not in theirsub]
    print('Number of mismatches with us:',len(mismatch))

    #Slice out corresponding columns.
    outtab = df.loc[keepSub,(['Subjects'] + cogtests)]

    #Replace names.
    newdict = {'CardSort_Unadj':'CardSort',
                'Flanker_Unadj':'Flanker',
                'ListSort_Unadj':'ListSort',
                'PicSeq_Unadj':'PicSeq',
                'PicVocab_Unadj':'PicVoc',
                'ProcSpeed_Unadj':'ProcSpeed',
                'ReadEng_Unadj':'ReadVoc',
                'PMAT24_A_CR':'PMAT24',
                'IWRD_TOT':'WordMem',
                'VSPLOT_TC':'LineOrient'}
    outtab.columns = pd.Series(outtab.columns).replace(newdict)

    #Output the table.
    subgroup = 'ourfa_subs'
    outpath = ('../outputs/c_cognition/'+subgroup+'/')
    os.makedirs(outpath,exist_ok=True)
    outfile = (outpath+'cogscores.csv')
    outtab.to_csv(outfile,index=False)
        