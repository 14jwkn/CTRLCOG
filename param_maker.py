#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Read HCP1200 subjects with resting-state fMRI scans based on the directory 
structure extracted from DataLad.
Output:
'r_allHCP_subjects.txt' Resting-state HCP subjects are listed at each line.

"""

import os
import numpy as np

#Extract all available subject files.
hcp_path = '../inputs/data/human-connectome-project-openaccess/HCP1200/'
final = sorted(os.listdir(hcp_path))
np.savetxt('r_allHCP_subjects.txt',final,fmt='%s',delimiter="\n")
    