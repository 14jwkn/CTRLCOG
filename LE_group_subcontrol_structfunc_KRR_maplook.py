# -*- coding: utf-8 -*-
"""
For a given k for clustering, type of SC normalization, and type and percentage of thresholding, 
find Dice coefficients between significant Haufe scores and PFIT/MD maps. Generate p-values
for Dice coefficients using permutation. Find Spearman's correlation to quantify the relationship
and R2 for quadratic regression to quantify the statistical relationship between Haufe scores and 
controllability/principal cortical gradient. Generate p-values for the R2.

Output:
'dc_'+cctrl+'_'+cstatetype+'_'+ccog+'_covha.csv' Dice coefficient between significant Haufe scores and PFIT/MD maps.
'dc_'+cctrl+'_'+cstatetype+'_'+ccog+'_covha_pval.csv' p-values for Dice coefficient between significant Haufe scores and PFIT/MD maps.
'dc_'+cctrl+'_'+cstatetype+'_'+ccog+'_covha_round.csv' Dice coefficient between significant Haufe scores and PFIT/MD maps rounded.
'dc_'+cctrl+'_'+cstatetype+'_'+ccog+'_covha_pthres.csv' Dice coefficient between significant Haufe scores and PFIT/MD maps rounded and thresholded for significance.
cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'.csv' Spearman correlation between Haufe scores and controllability/principal cortical gradient.
cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'_round.csv' Spearman correlation between Haufe scores and controllability/principal cortical gradient rounded.
cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'_'+outlist_lab+'.csv' R2 for quadratic regression between Haufe scores and controllability/principal cortical gradient.
cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'_'+outlist_lab+'_pval.csv' p-value for R2 for quadratic regression between Haufe scores and controllability/principal cortical gradient.
cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'_'+outlist_lab+'_round.csv' R2 for quadratic regression between Haufe scores and controllability/principal cortical gradient rounded. 
cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'_'+outlist_lab+'_pthres.csv' R2 for quadratic regression between Haufe scores and controllability/principal cortical gradient rounded and thresholded for significance. 

Usage: 
    LE_group_subcontrol_structfunc_KRR_maplook.py <k> <sctype> <threstype> <thresval> <controltype> <statetype> <septype> <nrep> <inner_k> <outer_k> <wantblock> <wantperm>
    
Arguments:
    
    <k> k for k-means clustering
    <sctype> SC type of normalization
    <threstype> Threshold type
    <thresval> Threshold value
    <controltype> What control types are used
    <statetype> What statetypes are used
    <septype> What cognitive separation types are used
    <nrep> Number of CV repetitions
    <inner_k> Inner K in K-fold CV hyperparameter search
    <outer_k> Outer K in K-fold CV
    <wantblock> Which analysis block to run
    <wantperm> What number of permuations to run

"""

import os, sys, time, random
import numpy as np
import pandas as pd
import nibabel as nib
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from enigmatoolbox.permutation_testing import spin_test
from enigmatoolbox.permutation_testing import quad_spin_test
from brainspace.utils.parcellation import map_to_labels
from brainspace.null_models import SpinPermutations
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from matplotlib import pyplot as plt
from docopt import docopt

#Permutation indices for regular permutation.
def perm_idx_create(nidx,nspecialp):

    #Set up seed and values.
    cogseed = 12345
    random.seed(cogseed)
    lorig = tuple(range(nidx))
    lnum = list(range(nidx))
    pset = set()
    pset.add(tuple(lnum))

    #Repeat adding until there are the desired number of non-repeated index sets.
    while len(pset) < (nspecialp+1):
        random.shuffle(lnum)
        pset.add(tuple(lnum))
    pset = list(pset)  
    pset.remove(lorig)
    return(pset)

#Dice calculator.
def dice_coef(ccomp,ccomp2):

  #Calculate DC based on ab = number of true matches, and a and b = number of true in each comparator.
  ab = (ccomp & ccomp2).sum()
  a = ccomp.sum()
  b = ccomp2.sum()
  dc = 2*(ab/(a+b))

  #Return dice coefficient.
  return dc

#Dice permutation.
def dice_perm(cdice,ccomp,ccomp2,pset):

    #Repeat and permute ccomp.
    pset2_idx = np.array(pset).T
    nspecialp = np.shape(pset2_idx)[1]
    pset2_val = np.zeros_like(pset2_idx)
    for coidx in range(nspecialp):
        pset2_val[:,coidx] = ccomp[pset2_idx[:,coidx]]
    
    #Conduct dice coefficient over each column.
    dc_permdist = np.apply_along_axis(dice_coef,0,pset2_val,ccomp2)   

    #Find the number greater than or equal to the dice value and divide by number
    #of permutations to get one-sided p-value. Add 1 to numerator and denominator
    #to place the true value in the permutation distribution as well.
    cdicep = pd.DataFrame(dc_permdist).ge(cdice,axis=1).sum(axis=0).add(1).div((nspecialp+1))[0]

    #Return dice p-value.
    return(cdicep)

#Spin test p-value alteration for two-sided.
def spin_alter(ccorr,cpdist):
    xy_pdist = cpdist[0:nspecialp]
    yx_pdist = cpdist[nspecialp:(nspecialp*2)]

    #Two-sided significance.
    xy_p = pd.DataFrame(xy_pdist).abs().ge(np.abs(ccorr),axis=1).sum(axis=0).add(1).div((nspecialp+1))[0]
    yx_p = pd.DataFrame(yx_pdist).abs().ge(np.abs(ccorr),axis=1).sum(axis=0).add(1).div((nspecialp+1))[0]

    #Average xy and yx p-values.
    alter_p = (xy_p+yx_p)/2

    return alter_p

#Generate spin permutations.
def spin_perm_gen(cvec,sphere_l,sphere_r,labels_l,labels_r,nspecialp):

    #Convert data to hemisphere and remove labels.
    nreg = len(cvec)
    nreg_h = int(nreg/2)
    cvec_r = cvec[:nreg_h].values
    cvec_l = cvec[nreg_h:nreg].values

    #Indexing of labels is 1 to 180 for right hemisphere and 181 to 360 
    #for left hemisphere. Turn 0s into NA. Convert to 0 to 179 for both.
    adj_labels_r = labels_r.copy().astype(float)
    adj_labels_r[adj_labels_r==0] = np.nan
    adj_labels_r = adj_labels_r - 1
    adj_labels_l = labels_l.copy().astype(float)
    adj_labels_l[adj_labels_l==0] = np.nan
    adj_labels_l = adj_labels_l - 181

    #Map parcellated data to vertex-level data for each hemisphere.
    data_r = map_to_labels(cvec_r,adj_labels_r,ignore_missing=True)
    data_l = map_to_labels(cvec_l,adj_labels_l,ignore_missing=True)
    
    #Fit the spin model.
    spins = SpinPermutations(n_rep=nspecialp,random_state=12345)
    spins.fit(sphere_l,sphere_r,parcellation=(labels_l,labels_r))

    #Spin the vector into n permutations.
    spin_mat = spins.randomize(cvec)

    #Return spin-permuted vectors.
    return spin_mat

#Quadratic regression.
def quad_reg(cx,cy):
   
    #Generate squared term and add to model.
    cx2 = cx**2
    full_x = np.column_stack([cx2,cx])

    #Fit model.
    cmodel = LinearRegression().fit(full_x,cy)

    #Extract coefficients and R2.
    b2 = cmodel.coef_[0]
    b1 = cmodel.coef_[1]
    b0 = cmodel.intercept_
    r2 = cmodel.score(full_x,cy)

    #Return.
    return (b2,b1,b0,r2)

#Two-sided test for an observation and a vector of permutations.
def two_perm(cobs,cperms):
    nperm = cperms.shape[0]
    cp = (np.sum(np.abs(cperms)>=np.abs(cobs))+1)/(nperm+1)
    return(cp)

#One-sided test for an observation and a vector of permutations.
def one_perm(cobs,cperms):
    nperm = cperms.shape[0]
    cp = (np.sum(cperms >= cobs)+1)/(nperm+1)
    return(cp)

#Reader for control and degree values.
def control_reader(nk,nroi,
                   basepath,scpath,sFCpath):

    #Read in the full Gramian error table to find states with errors.
    infile = (basepath+'/state_images/SC_sFC_dFC_gram.csv')
    gramall = pd.read_csv(infile,index_col=0,header=None)
    gramstates = gramall.index.values
    gramstates = gramstates[np.squeeze((gramall!=0).values)].tolist()
    ngram = len(gramstates)
        
    #Get FC AC, MC, and D.
    infile = (basepath+'ave_tab.csv')
    dFCave = pd.read_csv(infile,index_col=0,header=None)
    infile = (basepath+'mod_tab.csv')
    dFCmod = pd.read_csv(infile,index_col=0,header=None)
    infile = (basepath+'deg_tab.csv')
    dFCdeg = pd.read_csv(infile,index_col=0,header=None)

    #Generate FC predictor labels.
    predlabs = []
    for kidx in range(nk):
        for ridx in range(nroi):
            predlabs.append('s'+str(kidx+1)+'_r'+str(ridx+1))
    dFCave.columns = predlabs
    dFCmod.columns = predlabs
    dFCdeg.columns = predlabs

    #If the error list contains dFC states.
    for kidx in range(nk):
        errstate = ('s'+str(kidx+1))
        if (errstate in gramstates):
            
            #Read in the regions to remove.
            infile = (basepath+'/state_images/'+errstate+'_err.csv')
            errmat = pd.read_csv(infile,index_col=None,header=0)
            regrem = errmat.loc[:,'Region'].values.tolist()
            regrem = [(errstate+'_r'+str(x)) for x in regrem]
            newlabs = [x for x in predlabs if x not in regrem]

            #Remove those from the state.
            dFCave = dFCave.loc[:,newlabs]
            predlabs = dFCave.columns
            dFCmod = dFCmod.loc[:,newlabs]
            predlabs = dFCmod.columns
            dFCdeg = dFCdeg.loc[:,newlabs]
            predlabs = dFCdeg.columns

    #Read SC.
    infile = (scpath+'orig_ave_sc.csv')
    scave = pd.read_csv(infile,index_col=0,header=None)
    infile = (scpath+'orig_mod_sc.csv')
    scmod = pd.read_csv(infile,index_col=0,header=None)
    infile = (scpath+'deg_sc.csv')
    scdeg = pd.read_csv(infile,index_col=0,header=None)
    predlabs = [('sc_r'+str(ridx+1)) for ridx in range(nroi)]
    scave.columns = predlabs
    scmod.columns = predlabs
    scdeg.columns = predlabs

    #If the error list contains SC.
    errstate = 'SC'
    if (errstate in gramstates):
        
        #Read in the regions to remove.
        infile = (basepath+'/state_images/'+errstate+'_err.csv')
        errmat = pd.read_csv(infile,index_col=None,header=0)
        regrem = errmat.loc[:,'Region'].values.tolist()
        regrem = [('sc_r'+str(x)) for x in regrem]
        newlabs = [x for x in predlabs if x not in regrem]

        #Remove those from the state.
        scave = scave.loc[:,newlabs]
        scmod = scmod.loc[:,newlabs]
        scdeg = scdeg.loc[:,newlabs]

    #Read sFC.
    infile = (sFCpath+'ave_sFC.csv')
    sFCave = pd.read_csv(infile,index_col=0,header=None)
    infile = (sFCpath+'mod_sFC.csv')
    sFCmod = pd.read_csv(infile,index_col=0,header=None)
    infile = (sFCpath+'deg_sFC.csv')
    sFCdeg = pd.read_csv(infile,index_col=0,header=None)
    predlabs = [('sFC_r'+str(ridx+1)) for ridx in range(nroi)]
    sFCave.columns = predlabs
    sFCmod.columns = predlabs
    sFCdeg.columns = predlabs

    #If the error list contains sFC.
    errstate = 'sFC'
    if (errstate in gramstates):
        
        #Read in the regions to remove.
        infile = (basepath+'/state_images/'+errstate+'_err.csv')
        errmat = pd.read_csv(infile,index_col=None,header=0)
        regrem = errmat.loc[:,'Region'].values.tolist()
        regrem = [('sFC_r'+str(x)) for x in regrem]
        newlabs = [x for x in predlabs if x not in regrem]

        #Remove those from the state.
        sFCave = sFCave.loc[:,newlabs]
        sFCmod = sFCmod.loc[:,newlabs]
        sFCdeg = sFCdeg.loc[:,newlabs]
        sFCgr = sFCgr.loc[newlabs]
    
    #Split single matrices into states.
    ctrl_collect = {}
    outmats = [scave,scmod,scdeg,
              sFCave,sFCmod,sFCdeg]
    outkeys = ['ave.sc','mod.sc','deg.sc',
               'ave.sFC','mod.sFC','deg.sFC']
    nout = len(outkeys)
    for ouidx in range(nout):
        outkey = outkeys[ouidx]
        outmat = outmats[ouidx]
        ctrl_collect[outkey] = outmat
   
    #Split multi-matrices into states.
    outmats = [dFCave,dFCmod,dFCdeg]
    outctrls = ['ave','mod','deg']
    nout = len(outctrls)
    for cstate in [('s'+str(x+1)) for x in range(nk)]:
        for ouidx in range(nout):
            cctrl = outctrls[ouidx]
            inmat = outmats[ouidx]
            cpred = inmat.columns
            newpred = [x for x in cpred if cstate in x]
            outmat = inmat.loc[:,newpred]
            outkey = (cctrl+'.'+cstate)
            ctrl_collect[outkey] = outmat

    #Return each item.
    return (ctrl_collect)

#Reader for prediction accuracy values and p-values.
def predict_reader(basepath,nrep,inner_k,outer_k,
                   sctype,ctrltypes,nctrl,statetypes,nstatetype,septypes,nseptype,cog_septypes,
                   pstatetypes,npstatetype,pseptypes,npseptype,pcog_septypes):
    
    #Produce predictive values for each model.
    pred_collect = {}
    for ctidx in range(nctrl):
        cctrl = ctrltypes[ctidx]

        #For each state type.
        for stidx in range(nstatetype):
            cstatetype = statetypes[stidx]

            #Set split state list based on state type.
            csplits = []
            if 'SC' in cstatetype:
                csplits = csplits + ['sc']
            if 'sFC' in cstatetype:
                csplits = csplits + ['sFC']
            if 'dFC' in cstatetype:
                csplits = csplits + [('s'+str(x+1)) for x in range(nk)]

            #For each cognitive set.
            for csidx in range(nseptype):
                cseptype = septypes[csidx]

                #Set cognitive list based on cognitive set.
                ccoglist = cog_septypes[csidx]
                cncog = len(ccoglist)

                #For each cognitive variable.
                for cidx in range(cncog):
                    ccog = ccoglist[cidx]
                
                    #Read in the version and average across folds and repetitions.
                    infile = (basepath+'KRRXFS/'+cctrl+'_'+cstatetype+'_'+cseptype+
                                '_'+nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+
                                '/score_collect.h5')
                    inkey = ('/rep_covha_feat_'+ccog)
                    store = pd.HDFStore(infile,'r')
                    inmat = store.select(inkey)
                    store.close()
                    inmat = inmat.mean(axis=0)

                    #Split states.
                    for cstate in csplits:
                        cpred = inmat.index
                        newpred = [x for x in cpred if cstate in x]
                        outmat = inmat.loc[newpred]
                        outkey = (cctrl+'.'+cstatetype+'.'+ccog+'.'+cstate)
                        pred_collect[outkey] = outmat

    #Produce predictive values for each model.
    pval_collect = {}
    for ctidx in range(nctrl):
        cctrl = ctrltypes[ctidx]

        #For each state type.
        for stidx in range(npstatetype):
            cstatetype = pstatetypes[stidx]

            #Set split state list based on state type.
            csplits = []
            if 'SC' in cstatetype:
                csplits = csplits + ['sc']
            if 'sFC' in cstatetype:
                csplits = csplits + ['sFC']
            if 'dFC' in cstatetype:
                csplits = csplits + [('s'+str(x+1)) for x in range(nk)]

            #For each cognitive set.
            for csidx in range(npseptype):
                cseptype = pseptypes[csidx]

                #Set cognitive list based on cognitive set.
                ccoglist = pcog_septypes[csidx]
                cncog = len(ccoglist)

                #For each cognitive variable.
                for cidx in range(cncog):
                    ccog = ccoglist[cidx]
                
                    #Read in the version.
                    inpath = (basepath+'KRRXFS/controlcomp_statecomp_cogcomp'+
                                '_'+nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/')
                    infile = (inpath+cctrl+'_'+cstatetype+'_'+cseptype+'_TwoP_BHFDR_covha_feat.csv')
                    inmat = pd.read_csv(infile,index_col=0)
                    inmat = inmat.loc[:,ccog]
                    
                    #Split states.
                    for cstate in csplits:
                        cpred = inmat.index
                        newpred = [x for x in cpred if cstate in x]
                        outmat = inmat.loc[newpred]
                        outkey = (cctrl+'.'+cstatetype+'.'+ccog+'.'+cstate)
                        pval_collect[outkey] = outmat < 0.05
    
    #Return each item.
    return (pred_collect,pval_collect)

if __name__ == '__main__':
    __spec__ = None
    
    #Catches arguments.
    args = docopt(__doc__)
    k = args['<k>']
    sctype = args['<sctype>']
    threstype = args['<threstype>']
    thresval = args['<thresval>']
    controltype = args['<controltype>']
    statetype = args['<statetype>']
    septype = args['<septype>']
    nrep = args['<nrep>']
    inner_k = args['<inner_k>']
    outer_k = args['<outer_k>']
    wantblock = args['<wantblock>']
    wantperm = args['<wantperm>']
    print('Doing:',k,sctype,threstype,thresval,
          controltype,statetype,septype,
          nrep,inner_k,outer_k,wantblock,wantperm)

    #Set paths.
    subgroup = 'full'
    sc_subgroup = 'dr_full'
    basepath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                subgroup+'/'+k+'/SC_dFC/'+
                sc_subgroup+'/collect/'+threstype+'/'+thresval+'/')
    scpath = ('../outputs/d_SC/'+sc_subgroup+'/'+threstype+'/'+thresval+'/'+sctype+'/')
    sFCpath = ('../outputs/r_sFC/'+
                sc_subgroup+'/'+threstype+'/'+thresval+'/')
    explorepath = (basepath+'KRRXFS/'+controltype+'_'+statetype+'_'+septype+'_'+
                  nrep+'_'+inner_k+
                  '_'+outer_k+'_'+sctype+'/correxplore/specialp/')
    os.makedirs(explorepath,exist_ok=True)
    
    #Initialize basic values.
    nspecialp = int(wantperm)
    nroi = 360
    nk = int(k)
    with open('dr_full_intersect.txt') as f:
        sublist = [label.rstrip() for label in f] 
    sublist = list(map(int,sublist))
    nsub = len(sublist)
    corrtype_list = ['spearman','quadreg']
    ncorrtype = len(corrtype_list)
    full_statesplits = ['sc','sFC'] + [('s'+str(x+1)) for x in range(nk)]

    #Set types of interest to be compared.
    statetypes = ['SC_dFCcat']
    nstatetype = len(statetypes)
    ctrltypes = ['ave','mod']
    nctrl = len(ctrltypes)
    septypes = ['comCFAng']
    nseptype = len(septypes)
    cog_septypes = [['gCFA']]
    coglist = sum(cog_septypes,[])
    ncog = len(coglist)

    #Produce versions for p-values.
    pstatetypes = ['SC_dFCcat']
    npstatetype = len(pstatetypes)
    pseptypes = ['comCFAng']
    npseptype = len(pseptypes)
    pcog_septypes = [['gCFA','P24_CR','PV','gFngCFA','gCngCFA']]
    pcoglist = sum(pcog_septypes,[])
    pncog = len(pcoglist)

    #Read in spherical indices for spin test.
    inpath = '../inputs/data/hcp/HCP_S1200_GroupAvg_v1/'
    infile = (inpath+'S1200.L.sphere.32k_fs_LR.surf.gii')
    sphere_l = nib.load(infile).darrays[0].data
    infile = (inpath+'S1200.R.sphere.32k_fs_LR.surf.gii')
    sphere_r = nib.load(infile).darrays[0].data

    #Read in atlas annotations for spin test.
    # Add HCP-MMP1.0 atlas file: '../inputs/data/atlas/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii'
    # Do command: wb_command -cifti-separate Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii COLUMN -label CORTEX_RIGHT Q1-Q6_rh_labels.label.gii
    infile = '../inputs/data/atlas/Q1-Q6_lh_labels.label.gii'
    labels_l = nib.load(infile).darrays[0].data
    infile = '../inputs/data/atlas/Q1-Q6_rh_labels.label.gii'
    labels_r = nib.load(infile).darrays[0].data

    #Read in area thresholds of interest.
    #Generated PFIT regions by manually annotating each HCP-MMP1.0 atlas region based on the PFIT regions described by:
    #Genç, E., Metzen, D., Fraenz, C., Schlüter, C., Voelkle, M. C., Arning, L., ... & Kumsta, R. (2023). Structural architecture and brain network efficiency link polygenic scores to intelligence. Human brain mapping, 44(8), 3359-3376.
    #Generated MD regions by manually annotating each HCP-MMP1.0 atlas region based on the extended MD system described by:
    #Assem, M., Glasser, M. F., Van Essen, D. C., & Duncan, J. (2020). A domain-general cognitive core defined in multimodally parcellated human cortex. Cerebral Cortex, 30(8), 4361-4380.
    areathres_collect = []
    for cview in ['PFIT','MD']:
        inmat = pd.read_csv('../outputs/r_area/'+cview+'/thres_'+cview+'_glasser.csv')
        areathres_collect.append(inmat)
    areathres_collect = pd.concat(areathres_collect,axis=1)
    areathres_labs = areathres_collect.columns.to_list()
    nareathres = len(areathres_labs)

    #Read in control and degree values.
    ctrl_collect = control_reader(nk,nroi,
                                  basepath,scpath,sFCpath)
    
    #Read in predictive values and corresponding p-values.
    [pred_collect,pval_collect] = predict_reader(basepath,nrep,inner_k,outer_k,sctype,
                                                ctrltypes,nctrl,statetypes,nstatetype,septypes,nseptype,cog_septypes,
                                                pstatetypes,npstatetype,pseptypes,npseptype,pcog_septypes)

    #Read in gradient values and flip to desired sign.
    infile = ('../outputs/r_sFC/unclean/both/whole/dr_full/none/0/sFC_gradients.csv')
    sFCgr_all = pd.read_csv(infile,header=None)
    infile = ('../outputs/r_sFC/unclean/both/whole/dr_full/none/0/sFC_gradients_flip.csv')
    sFCgr_flip = pd.read_csv(infile,header=None).values.tolist()[0]
    ngr = len(sFCgr_flip)
    for gidx in range(ngr):
        if sFCgr_flip[gidx] == 'T':
            sFCgr_all.iloc[:,gidx] = -sFCgr_all.iloc[:,gidx]
    sFCgr_all.index = [('r'+str(ridx+1)) for ridx in range(nroi)]
    
    #Do various area DC with significant predictiveness.
    if wantblock == 'areaDC':
        print('Area DC.')
        outlab = 'areaDC'
        outpath = (explorepath+outlab+'/')
        os.makedirs(outpath,exist_ok=True)

        #Generate permutation indices.
        nidx = nroi
        pset = perm_idx_create(nidx,nspecialp)
    
        #Set up labels.
        cstatetype = 'SC_dFCcat'
        csplits = []
        if 'SC' in cstatetype:
            csplits = csplits + ['sc']
        if 'sFC' in cstatetype:
            csplits = csplits + ['sFC']
        if 'dFC' in cstatetype:
            csplits = csplits + [('s'+str(x+1)) for x in range(nk)]
        nsplit = len(csplits)
        septypeint = 'comCFAng'
        cogint = ['gCFA']
        ncogint = len(cogint)

        #Read.
        for cctrl in ctrltypes:
            print(cctrl)
            for ccog in cogint:
                print(ccog)

                #Collect DC and average.
                fulldc = pd.DataFrame(np.zeros((nsplit,nareathres)),index=csplits,columns=areathres_labs)
                fulldc_p = pd.DataFrame(np.zeros((nsplit,nareathres)),index=csplits,columns=areathres_labs)
                for cstate in csplits:
                    print(cstate)

                    #Get p-values and predictiveness.
                    inkey = (cctrl+'.'+cstatetype+'.'+ccog+'.'+cstate)
                    cpthres = pval_collect[inkey]
                    cpred = pred_collect[inkey]

                    #Require full 360 regions. Replace errors with zero for zero predictiveness.
                    cpred_full = pd.Series(np.zeros(nroi,dtype=bool),index=[(cstate+'_r'+str(x+1)) for x in range(nroi)])
                    cpred_full[cpred.index] = cpthres
                    cpthres = cpred_full.copy()
                    cpred_full = pd.Series(np.zeros(nroi),index=[(cstate+'_r'+str(x+1)) for x in range(nroi)])
                    cpred_full[cpred.index] = cpred
                    cpred = cpred_full.copy()

                    #Make copies and remove indices.
                    cpthres = cpthres.copy().reset_index(drop=True)
                    cpred = cpred.copy().reset_index(drop=True)
    
                    #For each area.
                    for carea in areathres_labs:

                        #Get area threshold, make copy, and remove indices.
                        ccomp2 = areathres_collect.loc[:,carea].copy().reset_index(drop=True)

                        #Find the dice coefficient with the total, and find permutation p.
                        ccomp = cpthres
                        cdice = dice_coef(ccomp,ccomp2)
                        fulldc.loc[cstate,carea] = cdice
                        cdicep = dice_perm(cdice,ccomp,ccomp2,pset)
                        fulldc_p.loc[cstate,carea] = cdicep

                #Round, threshold, and save.
                fulldc_r = round(fulldc,2)
                fulldc_t = fulldc_r.copy()
                fulldc_t[fulldc_p>=0.05] = np.nan
                outfile = (outpath+'dc_'+cctrl+'_'+cstatetype+'_'+ccog+'_covha.csv')
                fulldc.to_csv(outfile)
                outfile = (outpath+'dc_'+cctrl+'_'+cstatetype+'_'+ccog+'_covha_pval.csv')
                fulldc_p.to_csv(outfile)
                outfile = (outpath+'dc_'+cctrl+'_'+cstatetype+'_'+ccog+'_covha_round.csv')
                fulldc_r.to_csv(outfile)
                outfile = (outpath+'dc_'+cctrl+'_'+cstatetype+'_'+ccog+'_covha_pthres.csv')
                fulldc_t.to_csv(outfile)

    #Do group-average controllability relationship with predictiveness.
    if wantblock == 'ctrlpred':
        print('Control-predictiveness relationship.')
        outlab = 'ctrlpred'
        outpath = (explorepath+outlab+'/')
        os.makedirs(outpath,exist_ok=True)

        #Set up.
        cstatetype = 'SC_dFCcat'
        csplits = []
        if 'SC' in cstatetype:
            csplits = csplits + ['sc']
        if 'sFC' in cstatetype:
            csplits = csplits + ['sFC']
        if 'dFC' in cstatetype:
            csplits = csplits + [('s'+str(x+1)) for x in range(nk)]
        nsplit = len(csplits)
        septypeint = 'comCFAng'
        cogint = ['gCFA']
        ncogint = len(cogint)

        #Read.
        for corrtype in corrtype_list:
            for cctrl in ctrltypes:

                #Do comparisons.
                fullcorr = pd.DataFrame(np.zeros((nsplit,ncogint)),index=csplits,columns=cogint)
                fullcorr_r = pd.DataFrame(np.zeros((nsplit,ncogint)),index=csplits,columns=cogint)
                fullquad_r2 = pd.DataFrame(np.zeros((nsplit,ncogint)),index=csplits,columns=cogint)
                fullquad_r2_p = pd.DataFrame(np.zeros((nsplit,ncogint)),index=csplits,columns=cogint)
                pdats = []
                for ccog in cogint:

                    #Do each state.
                    for cstate in csplits:

                        #Index predictiveness.
                        inkey = (cctrl+'.'+cstatetype+'.'+ccog+'.'+cstate)
                        cpred = pred_collect[inkey]
                        
                        #Requires full 360 regions. For predictiveness region errors, replace with 0 for zero predictiveness.
                        cpred_full = pd.Series(np.zeros(nroi),index=[(cstate+'_r'+str(x+1)) for x in range(nroi)])
                        cpred_full[cpred.index] = cpred
                        
                        #Index comparator, generate group-average.
                        inkey = (cctrl+'.'+cstate)
                        cmat = ctrl_collect[inkey].mean()

                        #Requires full 360 regions. For controllability region errors, replace with 1 for what the error value was
                        #as this is what is generated for regions with Gramian errors.
                        cmat_full = pd.Series(np.ones(nroi),index=[(cstate+'_r'+str(x+1)) for x in range(nroi)])
                        cmat_full[cmat.index] = cmat
                        
                        #Combine.
                        fullmat = pd.concat((cpred_full,cmat_full),axis=1)
                        fullmat.columns = ['Predict','Comparison']
                        pdats.append(fullmat)

                        #If doing correlation.
                        if (corrtype == 'spearman'):

                            #Correlate.
                            ccorr = spearmanr(fullmat.iloc[:,0],fullmat.iloc[:,1]).statistic
                            ccorr_r = round(ccorr,2)
                            
                            #Add.  
                            fullcorr.loc[cstate,ccog] = ccorr
                            fullcorr_r.loc[cstate,ccog] = ccorr_r

                        #If doing quadratic regression.
                        elif (corrtype == 'quadreg'):

                            #Extract x and y.
                            cx = fullmat.iloc[:,0]
                            cy = fullmat.iloc[:,1]

                            #Do quadratic regression.
                            b2,b1,b0,r2 = quad_reg(cx,cy)

                            #Do spin test.
                            np.random.seed(12345)   
                            b2_p,b1_p,b0_p,r2_p = quad_spin_test(fullmat.iloc[:,0],fullmat.iloc[:,1],surface_name='fsa5',parcellation_name='glasser_360',
                                                    type=corrtype,n_rot=nspecialp,null_dist=True)

                            #Save.
                            fullquad_r2.loc[cstate,ccog] = r2
                            fullquad_r2_p.loc[cstate,ccog] = r2_p

                #Save correlation.
                if (corrtype == 'spearman'):

                    #Save values and rounded values.
                    outfile = (outpath+cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'.csv')
                    fullcorr.T.to_csv(outfile)
                    outfile = (outpath+cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'_round.csv')
                    fullcorr_r.T.to_csv(outfile)
                    
                #Save quadratic regression.
                elif (corrtype=='quadreg'):

                    #Extract.
                    outmat = fullquad_r2
                    outmat_p = fullquad_r2_p
                    outlist_lab = 'r2'

                    #Save matrix and p values.
                    outfile = (outpath+cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'_'+outlist_lab+'.csv')
                    outmat.T.to_csv(outfile)
                    outfile = (outpath+cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'_'+outlist_lab+'_pval.csv')
                    outmat_p.T.to_csv(outfile)

                    #Save rounded and thresholded values.
                    outfile = (outpath+cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'_'+outlist_lab+'_round.csv')
                    outmat_ro = round(outmat,2)
                    outmat_ro.T.to_csv(outfile)
                    outfile = (outpath+cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'_'+outlist_lab+'_pthres.csv')
                    outmat_th = outmat_ro.copy()
                    outmat_th[outmat_p >= 0.05] = np.nan
                    outmat_th.T.to_csv(outfile)

    #Do sFC gradient 1 relationship with predictiveness.
    if wantblock == 'grad1pred':
        print('sFC gradient 1-predictiveness relationship.')
        outlab = 'grad1pred'
        outpath = (explorepath+outlab+'/')
        os.makedirs(outpath,exist_ok=True)

        #Set up.
        cstatetype = 'SC_dFCcat'
        csplits = []
        if 'SC' in cstatetype:
            csplits = csplits + ['sc']
        if 'sFC' in cstatetype:
            csplits = csplits + ['sFC']
        if 'dFC' in cstatetype:
            csplits = csplits + [('s'+str(x+1)) for x in range(nk)]
        nsplit = len(csplits)
        septypeint = 'comCFAng'
        cogint = ['gCFA']
        ncogint = len(cogint)
        sFCgr = sFCgr_all.iloc[:,0]

        #Read.
        for corrtype in corrtype_list:
            for cctrl in ctrltypes:

                #Do comparisons.
                fullcorr = pd.DataFrame(np.zeros((nsplit,ncogint)),index=csplits,columns=cogint)
                fullcorr_r = pd.DataFrame(np.zeros((nsplit,ncogint)),index=csplits,columns=cogint)
                fullquad_r2 = pd.DataFrame(np.zeros((nsplit,ncogint)),index=csplits,columns=cogint)
                fullquad_r2_p = pd.DataFrame(np.zeros((nsplit,ncogint)),index=csplits,columns=cogint)
                pdats = []
                for ccog in cogint:

                    #Do each state.
                    for cstate in csplits:

                        #Index predictiveness.
                        inkey = (cctrl+'.'+cstatetype+'.'+ccog+'.'+cstate)
                        cpred = pred_collect[inkey]
                        
                        #Requires full 360 regions. For the predictiveness error, replace with 0 for zero predictiveness.
                        cpred_full = pd.Series(np.zeros(nroi),index=[(cstate+'_r'+str(x+1)) for x in range(nroi)])
                        cpred_full[cpred.index] = cpred
                        
                        #Index comparator. Get all original values.
                        cmat = sFCgr.copy()
                        cmat.index = [(cstate+'_'+x) for x in cmat.index]
                        cmat_full = cmat

                        #Combine.
                        fullmat = pd.concat((cpred_full,cmat_full),axis=1)
                        fullmat.columns = ['Predict','Comparison']
                        pdats.append(fullmat)

                        #If doing correlations.
                        if (corrtype == 'spearman'):

                            #Correlate.
                            ccorr = spearmanr(fullmat.iloc[:,0],fullmat.iloc[:,1]).statistic
                            ccorr_r = round(ccorr,2)
                            
                            #Add.
                            fullcorr.loc[cstate,ccog] = ccorr
                            fullcorr_r.loc[cstate,ccog] = ccorr_r
                        
                        #If doing quadratic regression.
                        elif (corrtype == 'quadreg'):

                            #Extract x and y.
                            cx = fullmat.iloc[:,0]
                            cy = fullmat.iloc[:,1]

                            #Do quadratic regression.
                            b2,b1,b0,r2 = quad_reg(cx,cy)

                            #Do spin test.
                            np.random.seed(12345)   
                            b2_p,b1_p,b0_p,r2_p = quad_spin_test(fullmat.iloc[:,0],fullmat.iloc[:,1],surface_name='fsa5',parcellation_name='glasser_360',
                                                    type=corrtype,n_rot=nspecialp,null_dist=True)

                            #Save.
                            fullquad_r2.loc[cstate,ccog] = r2
                            fullquad_r2_p.loc[cstate,ccog] = r2_p
                        
                #Save correlations.
                if (corrtype == 'spearman'):

                    #Save values and rounded values.
                    outfile = (outpath+cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'.csv')
                    fullcorr.T.to_csv(outfile)
                    outfile = (outpath+cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'_round.csv')
                    fullcorr_r.T.to_csv(outfile)
                
                #Save quadratic regression.
                elif (corrtype=='quadreg'):
                    
                    #Extract.
                    outmat = fullquad_r2
                    outmat_p = fullquad_r2_p
                    outlist_lab = 'r2'

                    #Save matrix and p values.
                    outfile = (outpath+cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'_'+outlist_lab+'.csv')
                    outmat.T.to_csv(outfile)
                    outfile = (outpath+cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'_'+outlist_lab+'_pval.csv')
                    outmat_p.T.to_csv(outfile)

                    #Save rounded and thresholded values.
                    outfile = (outpath+cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'_'+outlist_lab+'_round.csv')
                    outmat_ro = round(outmat,2)
                    outmat_ro.T.to_csv(outfile)
                    outfile = (outpath+cctrl+'_'+cstatetype+'_'+septypeint+'_'+corrtype+'_'+outlab+'_'+outlist_lab+'_pthres.csv')
                    outmat_th = outmat_ro.copy()
                    outmat_th[outmat_p >= 0.05] = np.nan
                    outmat_th.T.to_csv(outfile)
    print('Exploration done.')
