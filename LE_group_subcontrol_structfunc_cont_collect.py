# -*- coding: utf-8 -*-
"""
For a given k for clustering, type of SC normalization, type and percentage of 
thresholding, collect continuous time system controllability versions for each state
from all subjects.
Output:
cont_'+controltype+'_'+subject+'.csv' CSV file containing continuous controllability for all states.

Usage: 
    LE_group_subcontrol_structfunc_cont_collect.py <k> <sctype> <threstype> <thresval> <controltype> 
    
Arguments:

    <k> k for k-means clustering
    <sctype> SC type of normalization
    <threstype> Threshold type
    <thresval> Threshold value
    <controltype> Control type

"""

import os, sys, time, h5py
import numpy as np
import pandas as pd
from docopt import docopt

if __name__ == '__main__':
    __spec__ = None
    
    #Catches arguments.
    args = docopt(__doc__)
    k = args['<k>']
    sctype = args['<sctype>']
    threstype = args['<threstype>']
    thresval = args['<thresval>']
    controltype = args['<controltype>']
    print('Doing:',k,sctype,threstype,thresval,controltype,)
  
    #Set paths.
    subgroup = 'full'
    sc_subgroup = 'dr_full'
    dFCpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
               subgroup+'/'+k+'/SC_dFC/'+sc_subgroup+'/collect/'+threstype+'/'+thresval)
    os.makedirs(dFCpath,exist_ok=True)

    #Read subjects.
    if (sc_subgroup == 'dr_full'):
        subfile = 'dr_full_submain.txt'
    with open(subfile) as f:
        subjects = [subject.rstrip() for subject in f]
    nsub = len(subjects)
    
    #Initialize values.
    nroi = 360
    nk = int(k)
    klabs = ['s'+str(x+1) for x in range(nk)]
    statelabs = ['sc','sFC'] + klabs
    nstates = len(statelabs)

    #Generate labels for all states and regions.
    all_cols = []
    for cstate in statelabs:
        for croi in range(nroi):
            all_cols.append(cstate+'_r'+str(croi+1))
    nall_cols = len(all_cols)

    #Read in for each subject.
    co_ctrl = pd.DataFrame(np.zeros((nsub,nall_cols)),index=subjects,columns=all_cols)
    for csub in subjects:

        #Read all states.
        infile = (dFCpath+'/cont_'+controltype+'_'+csub+'.csv')
        inmat = pd.read_csv(infile)

        #For each state, divide it into the collector.
        for cstate in statelabs:
            cctrl = inmat.loc[:,cstate]
            cctrl_labs = [(cstate+'_r'+str(croi+1)) for croi in range(nroi)]
            cctrl.index = cctrl_labs
            co_ctrl.loc[csub,cctrl_labs] = cctrl

    #Save.
    outfile = (dFCpath+'/cont_'+controltype+'_all.csv')
    outmat = pd.DataFrame(co_ctrl,columns=statelabs)
    outmat.to_csv(outfile,header=True)
    print('Continuous done.')
