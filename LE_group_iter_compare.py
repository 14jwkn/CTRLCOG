#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find the best iteration of k-medians clustering for each k. Values with lower 
within-cluster sum of distance (WCS) scores are better.
Adapts code from:
https://github.com/adolphslab/HCP_MRI-behavior/blob/master/intelligence.ipynb
Output:
best_iter.csv For all ks, contains the best iteration.  

"""

import os, h5py
import numpy as np
import pandas as pd
from docopt import docopt

if __name__ == '__main__':
    __spec__ = None
    
    #Set output folder path. If output folder doesn't exist, creates it.
    subgroup = 'full'
    outpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/best_iter/'+subgroup+'/')
    os.makedirs(outpath,exist_ok=True)
   
    #Set output file names. 
    iterfile = outpath + 'best_iter.csv'
    
    #Go through each iteration to make labels.
    row_labs = []
    klist = ['2','3','4','5','6','7','8','9','10','11','12']
    for k in klist:
        row_labs.append(k)
                                
    #Create empty matrix to save best iterations for the group.
    best_mat = pd.DataFrame(np.zeros((1,len(row_labs))))
    best_mat.columns = row_labs
    best_mat.index = [subgroup]
    
    #Go through each iteration and discover the best ones.               
    for k in klist:  
        niters = 500
        curr_lab = k
        bestwcs = float('inf')
        bestitx = 0
        for itx in range(1,niters+1):
            inpath = ('../outputs/r_stateflex/statecalc_test/LE/'+
                      'ver_MATLAB/'+subgroup+'/'+k+'/')  
            infile = (inpath+'substat_'+str(itx)+'.h5')
            inkey = ('/'+subgroup)
            if not os.path.exists(infile):
                print(curr_lab,str(itx),'does not exist.')
                continue
            statstore = h5py.File(infile,'r')
            statmat = np.array(statstore[inkey]).T
            statstore.close()
            wcs = statmat[0,:].sum()
            if wcs < bestwcs:
                bestwcs = wcs
                bestitx = itx
        best_mat.loc[:,curr_lab] = bestitx

    #Save the table of best iterations.
    best_mat.to_csv(iterfile,header=True,index=True)
