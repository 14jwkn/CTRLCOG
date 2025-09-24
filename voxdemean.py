#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For each resting-state scan for each subject ID, subtract mean BOLD timeseries
voxelwise.
Output:
demean_rfMRI_REST*_*.dtseries.nii Demeaned fMRI files.

Usage: 
    r_voxdemean.py <subject> 
    
Arguments:

    <subject> Subject ID

"""
import os
import numpy as np
import nibabel as nb
from docopt import docopt
from nilearn import image

if __name__ == '__main__':
    __spec__ = None
    
    #Gets subject ID.
    arguments = docopt(__doc__)
    subject = arguments['<subject>']
    print('Doing:',subject)
    
    #For each run.
    for run in ['rfMRI_REST1_LR','rfMRI_REST1_RL','rfMRI_REST2_LR',
                'rfMRI_REST2_RL']:
        
        #Define paths. If the input file doesn't exist or the output file 
        #exists, continue to the next run.
        print('Doing:',run)
        inpath = '../inputs/data/fmri/'+subject+'/'
        infile = run+'_Atlas_MSMAll_hp2000_clean.dtseries.nii'
        outfile = 'demean_' + infile
        if (not os.path.exists(inpath+infile)) or os.path.exists(inpath+outfile):
            print('Infile does not exist or outfile exists.')
            continue
        
        #Read in the image as a numpy array.
        nimg = image.load_img(inpath+infile)
        
        #Convert it to a numpy array.
        npy = nimg.get_fdata()
        
        #Demean.
        new_npy = npy - np.mean(npy,axis=0)
        
        #Convert it to image.
        new_nimg = nb.Cifti2Image(new_npy,header=nimg.header,
                                 nifti_header=nimg.nifti_header)
        
        #Save the image.
        new_nimg.to_filename(inpath+outfile) 
        
