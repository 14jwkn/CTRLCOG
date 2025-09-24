#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read in all participants with ICA-FIX 3T resting-state data, then isolate 
participants who have all four scans. Use the HCP1200 data dictionary to
isolate participants who have all 10 cognitive test scores of interest. Retrieve
the mean relative FD for each participant for each run and exclude participants
with higher than 0.2 in any run. Find the intersection across subjects.
Output:
r_full_submain.txt Main list of subject IDs in each line.

"""

import os
import pandas as pd
import numpy as np

if __name__ == '__main__':
    __spec__ = None

    #Gets all participants with whole MSMAll-processed 3T REST1-2 mean timeseries.
    sublist = 'allHCP_subjects.txt'
    with open(sublist) as f:
            subjects = [subject.rstrip() for subject in f]
    filepath = '../outputs/r_meants/'
    rest1l = 'demean_rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv'
    rest1r = 'demean_rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv'
    rest2l = 'demean_rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv'
    rest2r = 'demean_rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv'
    files = [rest1l,rest1r,rest2l,rest2r]
    errorlist = []
    for sub in subjects:
        print('Doing:',sub)
        for file in files:
            subpath = filepath+sub+'/'
            if os.path.exists(subpath+file):
                meants = pd.read_csv(subpath+file,header=None)
                if (meants.shape[0] != 360) or (meants.shape[1] != 1200):
                    errorlist.append(sub)
            else:
                errorlist.append(sub)
    meants_errors = sorted(list(set(errorlist)))
    allrest = [x for x in subjects if x not in meants_errors]

    #Gets all participants with all cognitive scores of interest.
    HCPdict = pd.read_csv('../inputs/data/hcp/HCP1200_Data_Dictionary.csv',index_col=0)
    HCPabrv = HCPdict[['CardSort_Unadj',
                    'Flanker_Unadj',
                    'ListSort_Unadj',
                    'PicSeq_Unadj',
                    'PicVocab_Unadj',
                    'ProcSpeed_Unadj',
                    'ReadEng_Unadj',
                    'PMAT24_A_CR',
                    'IWRD_TOT',
                    'VSPLOT_TC']]
    cleanHCPabrv = HCPabrv.dropna(axis=0)
    allcog = cleanHCPabrv.index
    allcog = list(map(str,list(allcog)))

    #Get all the participants who have an average FD < 0.2 in every run.
    nsub = len(subjects)
    runlist = ['REST1_LR','REST1_RL','REST2_LR','REST2_RL']
    nrun = len(runlist)
    motmat = pd.DataFrame(np.zeros(nsub,nrun))
    motmat.index = subjects
    motmat.columns = runlist
    for csub in subjects:
        for run in runlist:
            infile = ('../inputs/data/motion/'+csub+'Movement_RelativeRMS.txt')
            with open(infile) as f:
                fd = [val.rstrip() for val in f]
            motmat.loc[sub,run] = np.mean(fd)
    motmat.replace('None',np.nan,inplace=True)          
    motmat = motmat.apply(pd.to_numeric)
    highmot = motmat.query('REST1_LR>=0.2 or REST1_RL>=0.2 or REST2_LR>=0.2 or REST2_RL>=0.2')
    lowmot = list(set(list(motmat.index)).difference(list(highmot.index)))
    lowmot.sort()
    allmot = list(map(str,lowmot))

    #Finds the intersection between the participants.
    complist = sorted(list(set.intersection(*map(set,[allcog,allrest,allmot]))))

    #Produce final list. 
    complist.sort()
    np.savetxt('r_full_submain.txt',complist,delimiter="\n",fmt="%s")    

