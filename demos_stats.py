# -*- coding: utf-8 -*-
"""
Get gender, race, age, years of education, handedness, g, gF, and gC and print
descriptive statistics for them. Save gender subject lists to do gender-restricted
analyses.
Output:
dr_full_intersect_msub.txt Male subject list.
dr_full_intersect_fsub.txt Female subject list.

Usage: 
    demos_stats.py 

"""

import os
import pandas as pd
import numpy as np

# Read original subject file.
subfile = 'dr_full_intersect.txt'
with open(subfile) as f:
    fullsub = [subject.rstrip() for subject in f]

# Read HCP dictionary, including open and restricted.
hcp_open = pd.read_csv('../inputs/data/hcp/HCP1200_Data_Dictionary.csv',index_col=0)
hcp_res = pd.read_csv('../inputs/data/hcp/RESTRICTED_HCP1200_Data_Dictionary.csv',index_col=0)
hcp_dict = pd.concat((hcp_open,hcp_res),axis=1)
hcp_dict.index = hcp_dict.index.astype(str)
hcp_dict = hcp_dict.loc[fullsub,:]

# Read cognitive scores.
infile = ('../outputs/c_cognition/full/pred_all.csv')
cog_scores = pd.read_csv(infile,index_col=0,header=0)
cog_scores.index = cog_scores.index.astype(str)
cog_scores = cog_scores.loc[fullsub,:]

# Divide into male and female subject list and save.
msub = hcp_dict.index[hcp_dict.loc[:,'Gender']=='M']
fsub = hcp_dict.index[hcp_dict.loc[:,'Gender']=='F']
np.savetxt('dr_full_intersect_msub.txt',msub,delimiter="\n",fmt="%s")
np.savetxt('dr_full_intersect_fsub.txt',fsub,delimiter="\n",fmt="%s")

# Find age, years of education, handedness mean, SD, min, and max.
cdem = 'Age_in_Yrs'
print(cdem,hcp_dict.loc[:,cdem].mean(),hcp_dict.loc[:,cdem].std(),
      hcp_dict.loc[:,cdem].min(),hcp_dict.loc[:,cdem].max())
cdem = 'SSAGA_Educ'
print(cdem,hcp_dict.loc[:,cdem].mean(),hcp_dict.loc[:,cdem].std(),
      hcp_dict.loc[:,cdem].min(),hcp_dict.loc[:,cdem].max())
cdem = 'Handedness'
print(cdem,hcp_dict.loc[:,cdem].mean(),hcp_dict.loc[:,cdem].std(),
      hcp_dict.loc[:,cdem].min(),hcp_dict.loc[:,cdem].max())

# Find gender and race proportions.
cdem = 'Gender'
print(cdem,hcp_dict.loc[:,cdem].value_counts(normalize=True))
print(cdem,hcp_dict.loc[:,cdem].value_counts())
cdem = 'Race'
print(cdem,hcp_dict.loc[:,cdem].value_counts(normalize=True))
print(cdem,hcp_dict.loc[:,cdem].value_counts())

# Find g, gF, and gC mean and SD.
ccog = 'gCFA'
print(ccog,cog_scores.loc[:,ccog].mean(),cog_scores.loc[:,ccog].std(),
      cog_scores.loc[:,ccog].min(),cog_scores.loc[:,ccog].max())
ccog = 'P24_CR'
print(ccog,cog_scores.loc[:,ccog].mean(),cog_scores.loc[:,ccog].std(),
      cog_scores.loc[:,ccog].min(),cog_scores.loc[:,ccog].max())
ccog = 'PV'
print(ccog,cog_scores.loc[:,ccog].mean(),cog_scores.loc[:,ccog].std(),
      cog_scores.loc[:,ccog].min(),cog_scores.loc[:,ccog].max())
