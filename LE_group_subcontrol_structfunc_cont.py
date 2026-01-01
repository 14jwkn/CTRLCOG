# -*- coding: utf-8 -*-
"""
For a given k for clustering, type of SC normalization, type and percentage of 
thresholding for the subject, calculate continuous time system controllability 
versions for each state.
Output:
cont_'+controltype+'_'+subject+'.csv' CSV file containing continuous controllability for all states.

Usage: 
    LE_group_subcontrol_structfunc_cont.py <k> <subject> <sctype> <threstype> <thresval> <controltype> 
    
Arguments:

    <k> k for k-means clustering
    <subject> Subject
    <sctype> SC type of normalization
    <threstype> Threshold type
    <thresval> Threshold value
    <controltype> Control type

"""

import os, sys, time, h5py
import numpy as np
import pandas as pd
from nctpy.utils import matrix_normalization
from nctpy.metrics import ave_control, modal_control
from docopt import docopt

if __name__ == '__main__':
    __spec__ = None
    
    #Catches arguments.
    args = docopt(__doc__)
    k = args['<k>']
    subject = args['<subject>']
    sctype = args['<sctype>']
    threstype = args['<threstype>']
    thresval = args['<thresval>']
    controltype = args['<controltype>']
    print('Doing:',k,subject,sctype,threstype,thresval,controltype,)
  
    #Set paths.
    subgroup = 'full'
    sc_subgroup = 'dr_full'
    dFCpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
               subgroup+'/'+k+'/SC_dFC/'+sc_subgroup+'/collect/'+threstype+'/'+thresval)
    os.makedirs(dFCpath,exist_ok=True)
    
    #Initialize values.
    nroi = 360
    nk = int(k)
    
    #For each state, read it in.
    klabs = ['s'+str(x+1) for x in range(nk)]
    statelabs = ['sc','sFC'] + klabs
    nstates = len(statelabs)
    statelist = []

    #SC.
    inpath = ('../outputs/d_SC/'+sc_subgroup+'/'+threstype+'/'+subject+'/'+
                thresval+'/'+sctype)
    infile = (inpath+'/sub_sc.h5')
    instore = h5py.File(infile,'r')
    inkey = ('/'+subject)
    inmat = np.array(instore[inkey]).T
    instore.close()
    statelist.append(inmat)

    #sFC.
    inpath = ('../outputs/r_sFC/'+sc_subgroup+'/'+subject+'/'+threstype+'/'+thresval)
    infile = (inpath+'/sub_sFC.h5')
    instore = h5py.File(infile,'r')
    inkey = ('/'+subject)
    inmat = np.array(instore[inkey]).T
    instore.close()
    statelist.append(inmat)

    #dFC.
    inpath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
            +subgroup+'/'+k+'/SC_dFC/'+sc_subgroup+'/'+subject+'/'+threstype+'/'+thresval)
    infile = (inpath+'/sub_dFC.h5')
    for kidx in range(nk):
        cstate = ('s'+str(kidx+1))
        instore = h5py.File(infile,'r')
        inkey = ('/'+subject+'_'+cstate)
        inmat = np.array(instore[inkey]).T
        instore.close()
        statelist.append(inmat)
  
    #Calculate the continuous versions.
    co_ctrl = np.zeros((nroi,nstates))
    for stidx in range(nstates):
        print(stidx+1)
        a_raw = statelist[stidx]
        csys = 'continuous'
        a_norm = matrix_normalization(a_raw,system=csys) 
        if (controltype=='ave'):
            cctrl = ave_control(a_norm,system=csys) 
            co_ctrl[:,stidx] = cctrl
    
    #Save.
    outfile = (dFCpath+'/cont_'+controltype+'_'+subject+'.csv')
    outmat = pd.DataFrame(co_ctrl,columns=statelabs)
    outmat.to_csv(outfile,header=True)
    print('Continuous done.')
