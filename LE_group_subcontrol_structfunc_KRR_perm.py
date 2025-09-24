# -*- coding: utf-8 -*-
"""
For a given k for clustering, type of SC normalization, and type and percentage of thresholding, fit models from CV with a single 
permutation using the same parameters from the original model. The permuted CV accuracy for the test set and Haufe scores
are both saved.
Output:
cperm_'+permidx+'.h5' h5 file containing the accuracy score and Haufe scores for a single permutation. '

Usage: 
    LE_group_subcontrol_structfunc_KRRXFS_perm.py <k> <sctype> <threstype> <thresval> <controltype> <statetype> <septype> <nrep> <inner_k> <outer_k> <nperm> <permidx>
    
Arguments:
    
    <k> k for k-means clustering
    <sctype> SC type of normalization
    <threstype> Threshold type
    <thresval> Threshold value
    <controltype> Control type
    <statetype> Which state 
    <septype> Cognitive variables batch
    <nrep> Number of CV repetitions
    <inner_k> Inner K in K-fold CV hyperparameter search
    <outer_k> Outer K in K-fold CV
    <nperm> Number of permutations
    <permidx> Current permutation ID

"""

import os, sys, time, h5py, random
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from himalaya.backend import set_backend
from himalaya.kernel_ridge import MultipleKernelRidgeCV
from himalaya.kernel_ridge import Kernelizer
from himalaya.kernel_ridge import ColumnKernelizer
from himalaya.kernel_ridge._kernels import pairwise_kernels
from himalaya.kernel_ridge._random_search import _decompose_kernel_ridge
from himalaya.scoring import r2_score, correlation_score
from himalaya.scoring import r2_score_split
from himalaya.scoring import r2_score_split_svd
from docopt import docopt

#Do cross-validation for given data_X, data_Y, and data_C.
def overall_crossval_perm(data_X,data_Y,data_C, #Data.
                          permvals, #Permutation values.
                          oldY,oldalph,oldgam,olddelt, #Old input Y values and fitted hyperparameters.
                          nrep,outer_k,outercv_test, #CV items
                          mykernel,slices,nspace, #Multi-KRR items
                          repseed, #Random items
                          xcon): #Preprocessing items

    #Set up copies to make sure.
    data_X = data_X.copy()
    data_Y = data_Y.copy()
    data_C = data_C.copy()
    ncog = np.shape(data_Y)[1]

    #Define all CV repetition seeds from the current one.
    np.random.seed(repseed)
    rrepcv_list = np.random.randint(1,12345,nrep).tolist()
    
    #Go through CV repetitions.
    nrk = nrep*outer_k
    nsp = ncog*nspace
    [nsub,nfeat] = data_X.shape
    rep_testacc = np.zeros((nrk,ncog))
    rep_ha_feat = np.zeros((ncog,nrk,nfeat))
    rep_covha_feat = np.zeros((ncog,nrk,nfeat))
    rep_ynormha_feat = np.zeros((ncog,nrk,nfeat))
    nrkidx = 0
    for ridx in range(nrep):
        
        #Set the seed for this CV repetition.
        start1 = time.time()
        crep_seed = rrepcv_list[ridx]
        print('Repetition:',(ridx+1))

        #Define outer CV seeds.
        np.random.seed(crep_seed)
        routcv_list = np.random.randint(1,12345,outer_k).tolist()

        #Extract outer CV indices for this repetition.
        outercollect = outercv_test[:,ridx]
 
        #Go through outer CV loops.
        for outidx in range(outer_k):

            #Extract for current repetition and fold and permute.
            crk = ('r'+str(ridx+1)+'_f'+str(outidx+1))
            data_Y_new = oldY[crk]
            data_Y_new = data_Y_new[permvals.astype(int),:]
            data_Y_new = pd.DataFrame(data_Y_new,index=data_Y.index,columns=data_Y.columns)
            alphorig = oldalph[crk][0]
            gamorig = oldgam[crk]
            deltorig = olddelt[crk]

            #Extract indices.
            train_index = (np.where(outercollect!=(outidx+1))[0]).tolist()
            test_index = (np.where(outercollect==(outidx+1))[0]).tolist()
            
            #Extract train and test, setting up copies to make sure.
            X_train, X_test = data_X.iloc[train_index,:], data_X.iloc[test_index,:]
            Y_train_new, Y_test_new = data_Y_new.iloc[train_index,:], data_Y_new.iloc[test_index,:]
            C_train, C_test = data_C.iloc[train_index,:], data_C.iloc[test_index,:]
            X_train = X_train.copy()
            X_test = X_test.copy()
            Y_train_new = Y_train_new.copy()
            Y_test_new = Y_test_new.copy()
            C_train = C_train.copy()
            C_test = C_test.copy()
            
            #Confound regression on X.
            if xcon:

                #For each variable in X train.
                nX = np.shape(X_train)[1]
                new_X_train = np.zeros((np.shape(X_train)))
                new_X_test = np.zeros((np.shape(X_test)))
                for vidx in range(nX):

                    #Extract.
                    c_X_train = X_train.iloc[:,vidx]
                    c_X_test = X_test.iloc[:,vidx]

                    #Train a model on training X and apply it to get residuals.
                    conlr = LinearRegression().fit(C_train,c_X_train)
                    pred_X_train = conlr.predict(C_train)
                    new_c_X_train = (c_X_train - pred_X_train)
                    new_X_train[:,vidx] = new_c_X_train
                    
                    #Apply the model to the test X too.
                    pred_X_test = conlr.predict(C_test)
                    new_c_X_test = (c_X_test - pred_X_test)
                    new_X_test[:,vidx] = new_c_X_test
                
                #Reformat.
                X_train_index = X_train.index
                X_train_columns = X_train.columns
                X_test_index = X_test.index
                X_test_columns = X_test.columns
                X_train = pd.DataFrame(new_X_train,index=X_train_index,columns=X_train_columns)
                X_test = pd.DataFrame(new_X_test,index=X_test_index,columns=X_test_columns)
            
            #Standardize X_test using X_train mean and std.
            X_train_index = X_train.index
            X_train_columns = X_train.columns
            X_test_index = X_test.index
            X_test_columns = X_test.columns
            scaler = StandardScaler()
            X_train = pd.DataFrame(scaler.fit_transform(X_train),
                                   index=X_train_index,columns=X_train_columns)
            X_test = pd.DataFrame(scaler.transform(X_test),
                                  index=X_test_index,columns=X_test_columns)
            
            #Convert to float32 as they do for faster computation.
            X_train_np = backend.asarray(X_train.values,dtype=backend.float32)
            X_test_np = backend.asarray(X_test.values,dtype=backend.float32)
            Y_train_np = backend.asarray(Y_train_new.values,dtype=backend.float32)
            Y_test_np = backend.asarray(Y_test_new.values,dtype=backend.float32)
            ntrain = np.shape(X_train_np)[0]

            #Create train-train and hybrid test-train kernels for each space.
            same_list = []
            hybrid_list = []
            for spidx in range(nspace):
                colidx = slices[spidx]
                X_train_sp = X_train_np[:,colidx]
                X_test_sp = X_test_np[:,colidx]
                csame = pairwise_kernels(X_train_sp,X_train_sp,metric=mykernel)
                same_list.append(csame)
                chybrid = pairwise_kernels(X_test_sp,X_train_sp,metric=mykernel)
                hybrid_list.append(chybrid)
            Ks_same = backend.stack(same_list)
            Ks_hybrid = backend.stack(hybrid_list)

            #Go through combinations of kernels weighted by gamma for each target. 
            dual_weights = np.zeros((ntrain,ncog),dtype=np.float32)
            for cogidx in range(ncog):
                cgamma = gamorig[:,cogidx]
                cK = (cgamma[:,None,None] * Ks_same).sum(0) #Same.
                
                #Calculate weights.
                used_alphas = np.array([alphorig[cogidx]],dtype=np.float32)
                decomlist = list(_decompose_kernel_ridge(cK,used_alphas,Ktest=None,
                                                    negative_eigenvalues="zeros",
                                                    n_alphas_batch=len(used_alphas),
                                                    method='eigh'))
                matrix = decomlist[0][0] 
                weights = backend.matmul(matrix,Y_train_np[:,cogidx,None])

                #Append weights for this target.
                dual_weights[:,cogidx] = weights[0,:,0]

            #Multiply weights by alphas.
            dual_weights *= alphorig

            #Use these train weights to extract test set predictions and accuracy. 
            chi = backend.matmul(Ks_hybrid,dual_weights)
            split_predictions = backend.exp(deltorig[:, None,:]) * chi
            pred_Y_test_new = split_predictions.sum(0)
            acc_val_new = r2_score(Y_test_np,pred_Y_test_new)
            rep_testacc[nrkidx,:] = acc_val_new

            #Use these train weights to extract train set predictions.
            chi = backend.matmul(Ks_same,dual_weights)
            split_predictions = backend.exp(deltorig[:, None,:]) * chi
            pred_Y_train_new = split_predictions.sum(0)

            #Derive Haufe scores again.
            test_Y_haufe = pd.DataFrame(pred_Y_train_new,
                                       index=Y_train_new.index,
                                       columns=Y_train_new.columns)
            X_haufe = X_train.copy()
            nY = np.shape(test_Y_haufe)[1]
            for vidx in range(nY):

                #Extract.
                Y_haufe = test_Y_haufe.iloc[:,vidx]

                #Do covariance, training X is already demeaned. 
                X_haufe_demean = X_haufe
                Y_haufe_demean = Y_haufe - Y_haufe.mean()
                XY_cov = np.dot(X_haufe_demean.T,Y_haufe_demean)/(np.shape(X_haufe_demean)[0]-1)
                rep_covha_feat[vidx,nrkidx,:] = XY_cov.T

            #Increment k-fold loop and repetition.
            outidx += 1
            nrkidx += 1
        
        #Display time.
        end1 = time.time()
        print('Repetition done:',end1-start1)

    #Display diagnostic outputs.
    over_testacc = np.mean(rep_testacc,axis=0)
    print('Accuracy:',over_testacc)
    
    #Return values.
    return (rep_testacc,
            rep_covha_feat)

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
    nperm = args['<nperm>']
    permidx = args['<permidx>']
    print('Doing:',k,sctype,threstype,
          thresval,controltype,statetype,septype,nrep,
          inner_k,outer_k,nperm,permidx)
    
    #Set backend.
    cbackend = 'numpy' #Can be numpy, torch, torch_cuda, cupy
    backend = set_backend(cbackend,on_error='warn')
    print('Backend:',cbackend)

    #Set paths.
    subgroup = 'full'
    sc_subgroup = 'dr_full'
    basepath = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
                subgroup+'/'+k+'/SC_dFC/'+
                sc_subgroup+'/collect/'+threstype+'/'+thresval+'/')
    scpath = ('../outputs/d_SC/'+sc_subgroup+'/'+threstype+'/'+thresval+'/'+sctype+'/')
    sFCpath = ('../outputs/r_sFC/'+sc_subgroup+'/'+threstype+'/'+thresval+'/')
    krrpath = (basepath+'KRRXFS/'+controltype+'_'+statetype+'_'+septype+'_'+
              nrep+'_'+inner_k+'_'+outer_k+'_'+sctype+'/')
    permpath = (krrpath+'/perm_'+nperm+'/permsep/')
    os.makedirs(permpath,exist_ok=True)
    
    #Initialize values.
    nroi = 360
    nk = int(k)
    inner_k = int(inner_k)
    outer_k = int(outer_k) 
    nrep = int(nrep)
    alpha_list = [0.00001,0.0001,0.001,0.004,0.007,0.01,0.04,0.07,0.1,0.4,0.7,1,1.5,2,2.5,3,3.5,4,5,10,15,20]
    nperm = int(nperm)
    xcon = True
    ycon = True
    xother_cont = ['Age']
    xother_cat_list = ['Gender']
    xother_cat = '|'.join(xother_cat_list)
    
    #Quit if outfile exists.
    outfile = (permpath+'cperm_'+permidx+'.h5')
    if os.path.exists(outfile):
        print('Files exist.')
        sys.exit()

    #Read in the full group Gramian error table to find states with errors.
    infile = (basepath+'/state_images/'+sctype+'_SC_sFC_dFC_gram.csv')
    gramall = pd.read_csv(infile,index_col=0,header=None)
    gramstates = gramall.index.values
    gramstates = gramstates[np.squeeze((gramall!=0).values)].tolist()
    ngram = len(gramstates)
        
    #Get FC predictor data.
    infile = (basepath+controltype+'_tab.csv')
    dFCall = pd.read_csv(infile,index_col=0,header=None)
    sublist = list(dFCall.index)
    nsub = len(sublist)

    #Generate FC predictor labels.
    predlabs = []
    for kidx in range(nk):
        for ridx in range(nroi):
            predlabs.append('s'+str(kidx+1)+'_r'+str(ridx+1))
    dFCall.columns = predlabs

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
            dFCall = dFCall.loc[:,newlabs]
            predlabs = dFCall.columns

    #Read SC.
    if not ('deg' in controltype):
        infile = (scpath+'orig_'+controltype+'_sc.csv')
    else:
        infile = (scpath+'deg_sc.csv')
    scmat = pd.read_csv(infile,index_col=0,header=None)
    predlabs = [('sc_r'+str(ridx+1)) for ridx in range(nroi)]
    scmat.columns = predlabs

    #If the error list contains SC.
    errstate = 'SC'
    if (errstate in gramstates):
        
        #Read in the regions to remove.
        infile = (basepath+'/state_images/'+sctype+'_'+errstate+'_err.csv')
        errmat = pd.read_csv(infile,index_col=None,header=0)
        regrem = errmat.loc[:,'Region'].values.tolist()
        regrem = [('sc_r'+str(x)) for x in regrem]
        newlabs = [x for x in predlabs if x not in regrem]

        #Remove those from the state.
        scmat = scmat.loc[:,newlabs]

    #Read sFC.
    infile = (sFCpath+controltype+'_sFC.csv')
    sFCmat = pd.read_csv(infile,index_col=0,header=None)
    predlabs = [('sFC_r'+str(ridx+1)) for ridx in range(nroi)]
    sFCmat.columns = predlabs

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
        sFCmat = sFCmat.loc[:,newlabs]
    
    #Extract sizes of each feature space.
    dFC_size = np.shape(dFCall)[1]
    sc_size = np.shape(scmat)[1]
    sFC_size = np.shape(sFCmat)[1]

    #If doing SC, sFC, and dFC.
    if statetype == 'SC_sFC_dFCcat':

        #Concatenate everything.
        predmat = pd.concat((scmat,sFCmat,dFCall),axis=1)

        #Kernelize each slice.
        feature_names = ['sc','sFC','dFC']
        n_features_list = [sc_size,sFC_size,dFC_size]
        start_and_end = np.concatenate([[0],np.cumsum(n_features_list)])
        slices = [slice(start,end) 
                  for start,end in zip(start_and_end[:-1],start_and_end[1:])]
        nspace = len(feature_names)
    
    #If doing SC and dFC.
    elif statetype == 'SC_dFCcat':

        #Concatenate everything.
        predmat = pd.concat((scmat,dFCall),axis=1)

        #Kernelize each slice.
        feature_names = ['sc','dFC']
        n_features_list = [sc_size,dFC_size]
        start_and_end = np.concatenate([[0],np.cumsum(n_features_list)])
        slices = [slice(start,end) 
                  for start,end in zip(start_and_end[:-1],start_and_end[1:])]
        nspace = len(feature_names)
    
    #If doing SC and sFC.
    elif statetype == 'SC_sFC':

        #Concatenate everything.
        predmat = pd.concat((scmat,sFCmat),axis=1)

        #Kernelize each slice.
        feature_names = ['sc','sFC']
        n_features_list = [sc_size,sFC_size]
        start_and_end = np.concatenate([[0],np.cumsum(n_features_list)])
        slices = [slice(start,end) 
                  for start,end in zip(start_and_end[:-1],start_and_end[1:])]
        nspace = len(feature_names)
    
    #If doing dFC.
    elif statetype == 'dFCcat':

        #Concatenate everything.
        predmat = dFCall

        #Kernelize each slice.
        feature_names = ['dFC']
        n_features_list = [dFC_size]
        start_and_end = np.concatenate([[0],np.cumsum(n_features_list)])
        slices = [slice(start,end) 
                  for start,end in zip(start_and_end[:-1],start_and_end[1:])]
        nspace = len(feature_names)
    
    #If doing SC.
    elif statetype == 'SC':

        #Extract state.
        predmat = scmat

        #Kernelize each slice.
        feature_names = ['sc']
        n_features_list = [sc_size]
        start_and_end = np.concatenate([[0],np.cumsum(n_features_list)])
        slices = [slice(start,end) 
                  for start,end in zip(start_and_end[:-1],start_and_end[1:])]
        nspace = len(feature_names)
    
    #If doing sFC.
    elif statetype == 'sFC':

        #Extract state.
        predmat = sFCmat

        #Kernelize each slice.
        feature_names = ['sFC']
        n_features_list = [sFC_size]
        start_and_end = np.concatenate([[0],np.cumsum(n_features_list)])
        slices = [slice(start,end) 
                  for start,end in zip(start_and_end[:-1],start_and_end[1:])]
        nspace = len(feature_names)

    #See if there are any NA, NaN, or Inf in data.
    if predmat.isnull().any().any() or np.isinf(predmat).any().any():
        print('Value errors.')
        outfile = (krrpath+'error.txt')
        outmat = pd.DataFrame(['NaN or Inf in data.'])
        outmat.to_csv(outfile,index=None,header=None)
        sys.exit()

    #Get cognitive and other predictor data.
    infile = ('../outputs/c_cognition/'+subgroup+'/pred_all.csv')
    othercog = pd.read_csv(infile,index_col=0,header=0)
    othercog = othercog.loc[sublist,:]

    #Extract cognitive and continuous confounds.
    if septype == 'comCFAng':
        coglist = ['gCFA','P24_CR','PV','gFngCFA','gCngCFA']
    elif septype == 'comPCAng':
        coglist = ['gPCA2','gFngPCA2','gCngPCA2']
    elif septype == 'comEFAng':
        coglist = ['gEFA','gFngEFA','gCngEFA']
    elif septype == 'comNT':
        coglist = ['NT','NF','NC','NFnNT','NCnNT']
    ncog = len(coglist)
    cogmat = othercog.loc[:,coglist]
    othercont = othercog.loc[:,xother_cont]

    #Read in dummy data for categorical confounds, then merge.
    infile = ('../outputs/c_cognition/'+subgroup+'/dum_sel.csv')
    dummat = pd.read_csv(infile,index_col=0,header=0)
    dummat = dummat.loc[sublist,:]
    othercat = dummat.filter(regex=xother_cat)
    othermat = pd.concat((othercont,othercat),axis=1)

    #Read in the CV fold indices.
    infile = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
              subgroup+'/6/SC_dFC/'+
              sc_subgroup+'/predict_splits/r'+
              str(nrep)+'_o'+str(outer_k)+'_i'+str(inner_k)+'_predict_splits.h5')
    
    #Read in outer CV test indices for each repetition.
    inkey = '/outercv_testidx'   
    store = h5py.File(infile,'r')
    inmat = np.array(store[inkey]).T
    store.close()
    outercv_test = inmat

    #Read in inner CV test indices for each outer CV training set for each repetition.
    innercv_test = []
    for ridx in range(nrep):
        inkey = ('/r'+str(ridx+1)+'_innercv_testidx')
        store = h5py.File(infile,'r')
        inmat = np.array(store[inkey]).T
        store.close()
        innercv_test.append(inmat)

    #Read in permutation indices.
    infile = ('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/'+
              subgroup+'/6/SC_dFC/'+
              sc_subgroup+'/predict_splits/'+str(nperm)+'_perm_splits.h5')
    inkey = ('/permidx')
    store = h5py.File(infile,'r')
    inmat = np.array(store[inkey]).T
    store.close()
    permvals = inmat
    permvals = permvals[:,(int(permidx)-1)]
    permvals = permvals - 1

    #Generate repetition and k list.
    rk_labs = []
    for ridx in range(nrep):
        for fidx in range(outer_k):
            rk_labs.append('r'+str(ridx+1)+'_f'+str(fidx+1))
    nrk = len(rk_labs)

    #Read in the alpha, gamma, delta, and collected Y values from the original models.
    infile = (krrpath+'score_collect.h5')
    instore = pd.HDFStore(infile,'r')
    oldalph = {}
    oldgam = {}
    olddelt = {}
    oldY = {}
    for nrkidx in range(nrk):
        crk = rk_labs[nrkidx]
        inkey = ('/alpha_'+crk)
        oldalph[crk] = instore.select(inkey).values
        inkey = ('/gamma_'+crk)
        oldgam[crk] = instore.select(inkey).values
        inkey = ('/delta_'+crk)
        olddelt[crk] = instore.select(inkey).values
        inkey = ('/y_'+crk)
        oldY[crk] = instore.select(inkey).values
    instore.close()

    #Set X, Y with permutation, and C matrices.
    data_X = predmat
    data_Y = cogmat
    data_C = othermat
    
    #Run pipeline.
    mykernel = 'cosine'
    repseed = 12345
    [rep_testacc,
     rep_ha_feat,rep_covha_feat,rep_ynormha_feat] = overall_crossval_perm(data_X,data_Y,data_C, #Data
                                                                          permvals, #Permutation items
                                                                          oldY,oldalph,oldgam,olddelt, #Old input Y values and fitted hyperparameters.
                                                                          nrep,outer_k,outercv_test, #CV items
                                                                          mykernel,slices,nspace, #Multi-KRR items
                                                                          repseed, #Random items
                                                                          xcon) #Preprocessing items
 
    #Generate cognition and space list.
    cogspace_labs = []
    for cogidx in range(ncog):
        for spaidx in range(nspace):
            cogspace_labs.append(coglist[cogidx]+'_'+feature_names[spaidx])
    
    #Generate feature labels.
    feat_labs = predmat.columns

    #Save overall.
    outstore = pd.HDFStore(outfile)
    outkey = ('/rep_testacc')
    outstore.put(outkey,pd.DataFrame(rep_testacc,index=rk_labs,columns=coglist),format='table')

    #Go through feature matrices.
    for cidx in range(ncog):
        
        #Extract.
        ccog = coglist[cidx]

        #Save version.
        outkey = ('/rep_covha_feat_'+ccog)
        outstore.put(outkey,
                     pd.DataFrame(rep_covha_feat[cidx,:,:],columns=feat_labs),
                     format='table')
    outstore.close()
    print('Permutation done.')
