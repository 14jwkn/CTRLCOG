#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the selected k, use the sorting matrix to translate
the original clustering to the universal indices across k.
Output:
uni_subcent.h5 Universal index clustering.

Usage: 
    LE_group_stateclass.py <k> 
    
Arguments:

    <k> k for k-clustering

"""

import os, h5py
import numpy as np
import pandas as pd
from docopt import docopt

if __name__ == '__main__':
    __spec__ = None
    
    #Catches command-line arguments.
    args = docopt(__doc__)
    k = args['<k>']
    print('Doing:',k)
    
    #Set up I/O. 
    subgroup = 'full'    
    outpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
               subgroup+'/'+k+'/') 
    os.makedirs(outpath,exist_ok=True)
    transfile = (outpath+'uni_subcent.h5')    

    #Read in the centroid translation table.
    classfile = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/visclass/'+
              subgroup+'/allk/'+subgroup+'_sortclass.csv')
    classmat = pd.read_csv(classfile,header=0,index_col=0) 
    classmat.index = [str(x) for x in classmat.index]
    classmat = classmat.loc[k,:] 
    classmat = classmat.dropna().astype('int')
    classmat.index = pd.to_numeric(classmat.index,errors='coerce')
    classmat = classmat.to_dict() 
    classmat = {v: k for k, v in classmat.items()}     
      
    #Read in the best clustering.  
    bestfile = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/'+
              'group/best_iter/'+subgroup+'/best_iter.csv')
    bestiter = pd.read_csv(bestfile,header=0,index_col=0)
    bestlabel = k
    iteration = str(bestiter.loc[subgroup,bestlabel])
    infile = (outpath+'subclust_'+iteration+'.h5')
    inkey = ('/'+subgroup)
    instore = h5py.File(infile,'r')
    clustmat = np.array(instore[inkey],dtype='int').T
    instore.close()
    clustmat = pd.Series(np.squeeze(clustmat))
    
    #Translate.
    transmat = clustmat.map(classmat)
        
    #Save translation table.
    outstore = pd.HDFStore(transfile)
    outkey = ('/'+subgroup)
    outstore.put(outkey,transmat,format='table')
    outstore.close()  
    print('Saved.')   
      