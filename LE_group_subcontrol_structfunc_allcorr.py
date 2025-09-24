# -*- coding: utf-8 -*-
"""
Find the correlations between g, gF, gC, gF-g, and gC-g for the main estimate 
with CFA. Also find the correlations for estimates based on EFA, PCA, and the NIH toolbox.
Output:
corrmethod+'_'+cset+'_cogcorr.csv' Correlations between a set of cognitive variables.
corrmethod+'_'+cset+'_cogcorr_round.csv' Correlations rounded.

"""

import os, sys, h5py
import numpy as np
import pandas as pd

if __name__ == '__main__':
    __spec__ = None
    
    #Set outpath.
    subgroup = 'full'
    sc_subgroup = 'dr_full'
    cogpath = ('../outputs/c_cognition/'+sc_subgroup+'/')
    
    #Set parameters.
    corrmethod = 'pearson'
    
    #Set subject list.
    subfile = 'dr_full_intersect.txt'
    with open(subfile) as f:
        sublist = [subject.rstrip() for subject in f]
    
    #Read cognitive data.
    infile = ('../outputs/c_cognition/'+subgroup+'/pred_all.csv')
    othercog = pd.read_csv(infile,index_col=0,header=0)
    othercog.index = othercog.index.astype(str)

    #Cognitive correlations for main five cognitive variables and supplementary.
    cogsets = ['comCFAng','comEFAng','comPCAng','comNTng','comG']
    ncogsets = len(cogsets)
    cogset_cog = [['gCFA','P24_CR','PV','gFngCFA','gCngCFA'],
                  ['gEFA','P24_CR','PV','gFngEFA','gCngEFA'],
                  ['gPCA','P24_CR','PV','gFngPCA','gCngPCA'],
                  ['NT','NF','NC','NFnNT','NCnNT'],
                  ['gCFA','gEFA','gPCA','NT']]
    for coidx in range(ncogsets):
        cset = cogsets[coidx]
        coglist = cogset_cog[coidx]
        cogmat = othercog.loc[sublist,coglist]
        cogcorr = cogmat.corr(method=corrmethod)
        cogcorr_round = round(cogcorr,2)
        outfile = (cogpath+corrmethod+'_'+cset+'_cogcorr.csv')
        cogcorr.to_csv(outfile)  
        outfile = (cogpath+corrmethod+'_'+cset+'_cogcorr_round.csv')
        cogcorr_round.to_csv(outfile)  
