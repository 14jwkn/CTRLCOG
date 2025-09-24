#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Do PCA to extract the first PC score and loadings to use as gPCA. 
Output:
PCA_load.csv PCA loadings.
PCA_scores.csv PCA scores.

"""

import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
    
if __name__ == '__main__':
    __spec__ = None

    #Set subjects.
    subfile = 'r_full_submain.txt'
    with open(subfile) as f:
        oursub = [subject.rstrip() for subject in f]
    
    #Retrieve cognitive scores.
    openpath = '../inputs/data/hcp/HCP1200_Behavioral.csv'
    df = pd.read_csv(openpath,header=0)
    wantcog = ['CardSort_Unadj',
            'Flanker_Unadj',
            'ListSort_Unadj',
            'PicSeq_Unadj',
            'PicVocab_Unadj',
            'ProcSpeed_Unadj',
            'ReadEng_Unadj',
            'PMAT24_A_CR',
            'IWRD_TOT',
            'VSPLOT_TC']
    df = df.loc[df.loc[:,'Subject'].apply(str).isin(oursub),['Subject']+wantcog]

    #Slice out corresponding columns and convert it to numpy for processing.
    ngtable = df.loc[:,wantcog].values

    #Z-score standardization of scores.
    mgtable = ngtable - np.nanmean(ngtable,axis=0,keepdims=True)
    zgtable = mgtable / np.nanstd(ngtable,axis=0,ddof=1,keepdims=True)
    
    #Conduct PCA to get principal component scores.
    pca = PCA()
    pcscore = pca.fit_transform(zgtable)

    #Verifies PCA validity, finds factor loadings and explained variance of 
    #each component. Sets the names of the components.
    components = ['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']
    
    #Computes factor loadings and creates data frame from it.
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    loading_mat = pd.DataFrame(loadings,columns=components,index=wantcog)
   
    #Given that all the loadings for each cognitive test are negative for PC1
    #and a higher score means better performance for each test, we flip the 
    #sign on the PC1 loadings for better interpretability - a higher score 
    #means better performance.
    pcscore[:,0] = -pcscore[:,0]
    outtab = pd.DataFrame(pcscore,columns=components,index=oursub)
    outtab.index.name = 'Subject'

    #Output the table.
    subgroup = 'full'
    outpath = ('../outputs/c_cognition/'+subgroup+'/')
    os.makedirs(outpath,exist_ok=True)
    loading_mat.to_csv(outpath+'PCA_load.csv',header=True,index=True) 
    outtab.to_csv(outpath+'PCA_scores.csv',header=True,index=True)
        