# -*- coding: utf-8 -*-
"""
For a given k for clustering, type of SC normalization, and type and percentage of thresholding, fit a model with CV.
This model predicts each cognitive variable in the given batch using either average controllability, modal controllability,
or degree calculated from the set of state matrices specified. CV is done for the given repetitions, and inner and outer
fold number. To facilitate permutation testing with the same parameters set, for each fold, the Y variable with transformations
applied prior to model fitting, alpha, beta, and gamma are all saved. Next, the CV accuracy for the test set, the CV accuracy
for the inner CV for the best hyperparameter set, the decomposed CV accuracy for the test set for each feature band, and the
decomposed CV accuracy after demeaning of the test set according to the decomposition assumptions are all saved. Lastly, the
Haufe scores are also saved.
Output:
score_collect.h5 h5 file containing parameters to conduct permutation testing, accuracy scores, and Haufe scores.

Usage: 
    LE_group_subcontrol_structfunc_KRR_score.py <k> <sctype> <threstype> <thresval> <controltype> <statetype> <septype> <nrep> <inner_k> <outer_k>
    
Arguments:
    
    <k> k for k-means clustering
    <sctype> SC type of normalization
    <threstype> Threshold type
    <thresval> Threshold value
    <controltype> Control type
    <statetype> Which states are included in the model
    <septype> Cognitive variables batch
    <nrep> Number of CV repetitions
    <inner_k> Inner K in K-fold CV hyperparameter search
    <outer_k> Outer K in K-fold CV

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
from docopt import docopt

#Do cross-validation for given data_X, data_Y, and data_C.
def overall_crossval(data_X,data_Y,data_C, #Data
                    nrep,outer_k,inner_k,outercv_test,innercv_test, #CV items
                    alpha_list,mykernel,feature_names,slices,nspace, #Multi-KRR items
                    repseed, #Random items
                    xcon,ycon): #Preprocessing items

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
    rep_collect_Y = np.zeros((nrk,nsub,ncog))
    rep_alpha = np.zeros((nrk,ncog),dtype=backend.float32)
    rep_gamma = np.zeros((nrk,nspace,ncog),dtype=backend.float32)
    rep_delta = np.zeros((nrk,nspace,ncog),dtype=backend.float32)
    rep_testacc = np.zeros((nrk,ncog))
    rep_inneracc = np.zeros((nrk,ncog))
    rep_splitacc = np.zeros((nrk,nsp))
    rep_splitacc_strict = np.zeros((nrk,nsp))
    rep_covha_feat = np.zeros((ncog,nrk,nfeat))
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

            #Extract indices.
            train_index = (np.where(outercollect!=(outidx+1))[0]).tolist()
            test_index = (np.where(outercollect==(outidx+1))[0]).tolist()
            
            #Extract train and test, setting up copies to make sure.
            X_train, X_test = data_X.iloc[train_index,:], data_X.iloc[test_index,:]
            Y_train, Y_test = data_Y.iloc[train_index], data_Y.iloc[test_index]
            C_train, C_test = data_C.iloc[train_index,:], data_C.iloc[test_index,:]
            X_train = X_train.copy()
            X_test = X_test.copy()
            Y_train = Y_train.copy()
            Y_test = Y_test.copy()
            C_train = C_train.copy()
            C_test = C_test.copy()

            #Confound regression on Y.
            if ycon:

                #For each variable in Y train.
                nY = np.shape(Y_train)[1]
                new_Y_train = np.zeros((np.shape(Y_train)))
                new_Y_test = np.zeros((np.shape(Y_test)))
                for vidx in range(nY):

                    #Extract.
                    c_Y_train = Y_train.iloc[:,vidx]
                    c_Y_test = Y_test.iloc[:,vidx]

                    #Train a model on training X and apply it to get residuals.
                    conlr = LinearRegression().fit(C_train,c_Y_train)
                    pred_Y_train = conlr.predict(C_train)
                    new_c_Y_train = (c_Y_train - pred_Y_train)
                    new_Y_train[:,vidx] = new_c_Y_train
                    
                    #Apply the model to the test X too.
                    pred_Y_test = conlr.predict(C_test)
                    new_c_Y_test = (c_Y_test - pred_Y_test)
                    new_Y_test[:,vidx] = new_c_Y_test

                #Reformat.
                Y_train_index = Y_train.index
                Y_train_columns = Y_train.columns
                Y_test_index = Y_test.index
                Y_test_columns = Y_test.columns
                Y_train = pd.DataFrame(new_Y_train,index=Y_train_index,columns=Y_train_columns)
                Y_test = pd.DataFrame(new_Y_test,index=Y_test_index,columns=Y_test_columns)
            
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
            
            #Save train-test full altered Y for permutation purposes.
            collect_Y = np.zeros((nsub,ncog))
            collect_Y[train_index,:] = Y_train.copy()
            collect_Y[test_index,:] = Y_test.copy()
            rep_collect_Y[nrkidx,:,:] = collect_Y
            
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
            
            #Kernelize each slice in pipeline.
            preprocess_pipeline = make_pipeline(Kernelizer(kernel=mykernel))
            kernelizers_tuples = [(name,preprocess_pipeline,slice_)
                                 for name,slice_ in zip(feature_names,slices)]
            column_kernelizer = ColumnKernelizer(kernelizers_tuples)

            #Extract inner CV indices for this outer CV fold.
            innercollect = (innercv_test[ridx])[:,outidx]
            innercollect = innercollect[(innercollect)!=0]

            #Generate inner CV indices.
            inner_kf = []
            for inidx in range(inner_k):

                 #Extract indices.
                 inner_train_index = (np.where(innercollect!=(inidx+1))[0])
                 inner_test_index = (np.where(innercollect==(inidx+1))[0])
                 inner_kf.append((inner_train_index,inner_test_index))
            
            #Initiate solver with seed for this outer CV fold.
            solver = 'random_search'
            outer_seed = routcv_list[outidx]
            solver_params = dict(alphas=np.array(alpha_list),
                                 progress_bar=False,
                                 score_func=r2_score)
            cmodel = MultipleKernelRidgeCV(kernels='precomputed',
                                           solver=solver,solver_params=solver_params,
                                           cv=inner_kf,random_state=outer_seed)

            #Do inner CV to get best hyperparameters to fit outer fold.
            pipeline = make_pipeline(column_kernelizer,cmodel)
            fitted = pipeline.fit(X_train,Y_train)

            #Extract original alpha, delta, and gamma.
            alphorig = fitted[1].best_alphas_
            gamorig = fitted[1].best_gammas_
            deltorig = fitted[1].deltas_
            rep_alpha[nrkidx,:] = alphorig
            rep_gamma[nrkidx,:,:] = gamorig
            rep_delta[nrkidx,:,:] = deltorig

            #Check whether lambda limit was reached.
            if (np.max(alphorig) == max(alpha_list)) or (np.min(alphorig) == min(alpha_list)):
                print('Lambda limit reached:',alphorig)

            #Get test accuracy value for this fold. Uses COD metric, higher is better.
            acc_val = backend.to_numpy(fitted.score(X_test,Y_test))
            rep_testacc[nrkidx,:] = acc_val

            #Get the CV accuracy score in the inner CV for the best hyperparameter set.
            #Uses COD metric, higher is better.
            acc_inner = np.max(backend.to_numpy(fitted[1].cv_scores_),axis=0)
            rep_inneracc[nrkidx,:] = acc_inner

            #Split this accuracy value into component parts using product measure.
            pred_Y_test = fitted.predict(X_test,split=True)
            split_scores = backend.to_numpy(r2_score_split(Y_test,pred_Y_test))
            split_format = pd.DataFrame(split_scores).T.stack().reset_index().iloc[:,2]
            rep_splitacc[nrkidx,:] = split_format

            #Split this accuracy value into component parts using product measure.
            #after demeaning Y test which is part of its assumptions.
            Y_test_demean = Y_test.copy()
            Y_test_demean -= Y_test.mean(0)
            split_scores = backend.to_numpy(r2_score_split(Y_test_demean,pred_Y_test))
            split_format = pd.DataFrame(split_scores).T.stack().reset_index().iloc[:,2]
            rep_splitacc_strict[nrkidx,:] = split_format

            #Fit on X_train and calculate from X_train for feature scores.
            pred_Y_train = pd.DataFrame(fitted.predict(X_train))
            X_haufe = X_train.copy()
            all_Y_haufe = pred_Y_train.copy()

            #Go through each cognitive variable to generate Haufe scores.
            nY = np.shape(all_Y_haufe)[1]
            for vidx in range(nY):

                #Extract.
                Y_haufe = all_Y_haufe.iloc[:,vidx]

                #Do covariance. Training X is already demeaned. 
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
    minlam = rep_alpha.min().min()
    maxlam = rep_alpha.max().max()
    lamlab = (str(minlam)+'-'+str(maxlam))
    print('Lambda:',lamlab)
    print('Accuracy:',over_testacc)
    
    #Return values.
    return (rep_collect_Y,rep_alpha,rep_gamma,rep_delta,
            rep_testacc,rep_inneracc,
            rep_splitacc,rep_splitacc_strict,
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
    print('Doing:',k,sctype,threstype,thresval,
          controltype,statetype,septype,
          nrep,inner_k,outer_k)
    
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
    os.makedirs(krrpath,exist_ok=True)
    
    #Initialize values.
    nroi = 360
    nk = int(k)
    inner_k = int(inner_k)
    outer_k = int(outer_k) 
    nrep = int(nrep)
    alpha_list = [0.00001,0.0001,0.001,0.004,0.007,0.01,0.04,0.07,0.1,0.4,0.7,1,1.5,2,2.5,3,3.5,4,5,10,15,20]
    xcon = True
    ycon = True
    xother_cont = ['Age']
    xother_cat_list = ['Gender']
    xother_cat = '|'.join(xother_cat_list)
    
    #Outfile define.
    outfile = (krrpath+'score_collect.h5')

    #Read in the full group Gramian error table to find states with errors.
    infile = (basepath+'/state_images/'+sctype+'_SC_sFC_dFC_gram.csv')
    gramall = pd.read_csv(infile,index_col=0,header=None)
    gramstates = gramall.index.values
    gramstates = gramstates[np.squeeze((gramall!=0).values)].tolist()
    ngram = len(gramstates)
        
    #Get dFC predictor data.
    infile = (basepath+controltype+'_tab.csv')
    dFCall = pd.read_csv(infile,index_col=0,header=None)
    sublist = list(dFCall.index)
    nsub = len(sublist)

    #Generate dFC predictor labels.
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

    #Extract cognitive and continuous confounds for the batch.
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
              subgroup+'/6/SC_dFC/'+sc_subgroup+'/predict_splits/r'+
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
        inkey = ('r'+str(ridx+1)+'_innercv_testidx')
        store = h5py.File(infile,'r')
        inmat = np.array(store[inkey]).T
        store.close()
        innercv_test.append(inmat)
    
    #Set X, Y, and C matrices.
    data_X = predmat
    data_Y = cogmat
    data_C = othermat

    #Run pipeline.
    mykernel = 'cosine'
    repseed = 12345
    [rep_collect_Y,rep_alpha,rep_gamma,rep_delta,
    rep_testacc,rep_inneracc,
    rep_splitacc,rep_splitacc_strict,
    rep_covha_feat] = overall_crossval(data_X,data_Y,data_C, #Data
                                        nrep,outer_k,inner_k,outercv_test,innercv_test, #CV items
                                        alpha_list,mykernel,feature_names,slices,nspace, #Multi-KRR items
                                        repseed, #Random items
                                        xcon,ycon) #Preprocessing items
    
    #Generate repetition and k list.
    rk_labs = []
    for ridx in range(nrep):
        for fidx in range(outer_k):
            rk_labs.append('r'+str(ridx+1)+'_f'+str(fidx+1))
    nrk = len(rk_labs)

    #Generate cognition and space list.
    cogspace_labs = []
    for cogidx in range(ncog):
        for spaidx in range(nspace):
            cogspace_labs.append(coglist[cogidx]+'_'+feature_names[spaidx])
    
    #Generate feature labels.
    feat_labs = predmat.columns

    #Initialize.
    outstore = pd.HDFStore(outfile)

    #Save parameters for each fold.
    for nrkidx in range(nrk):
        crk = rk_labs[nrkidx]

        #Save Y, alpha, gamma, and delta.
        cout = pd.DataFrame(rep_collect_Y[nrkidx,:,:],index=sublist,columns=coglist)
        outkey = ('/y_'+crk)
        outstore.put(outkey,cout,format='table')
        cout = pd.DataFrame(rep_alpha[nrkidx,:,None].T,columns=coglist)
        outkey = ('/alpha_'+crk)
        outstore.put(outkey,cout,format='table')
        cout = pd.DataFrame(rep_gamma[nrkidx,:,:],index=feature_names,columns=coglist)
        outkey = ('/gamma_'+crk)
        outstore.put(outkey,cout,format='table')
        cout = pd.DataFrame(rep_delta[nrkidx,:,:],index=feature_names,columns=coglist)
        outkey = ('/delta_'+crk)
        outstore.put(outkey,cout,format='table')

    #Save overall.
    outkey = ('/rep_testacc')
    outstore.put(outkey,pd.DataFrame(rep_testacc,index=rk_labs,columns=coglist),format='table')
    outkey = ('/rep_inneracc')
    outstore.put(outkey,pd.DataFrame(rep_inneracc,index=rk_labs,columns=coglist),format='table')
    outkey = ('/rep_splitacc')
    outstore.put(outkey,pd.DataFrame(rep_splitacc,index=rk_labs,columns=cogspace_labs),format='table')
    outkey = ('/rep_splitacc_strict')
    outstore.put(outkey,pd.DataFrame(rep_splitacc_strict,index=rk_labs,columns=cogspace_labs),format='table')

    #Go through feature matrices.
    for cidx in range(ncog):
        
        #Extract.
        ccog = coglist[cidx]

        #Save version.
        if int(k) <= 10:
            outkey = ('/rep_covha_feat_'+ccog)
            outstore.put(outkey,
                        pd.DataFrame(rep_covha_feat[cidx,:,:],columns=feat_labs),
                        format='table')
        else:
            outkey = ('/rep_covha_feat_'+ccog)
            outstore.put(outkey,
                        pd.DataFrame(rep_covha_feat[cidx,:,:],columns=feat_labs).T,
                        format='table')
    outstore.close()
    print('Scoring done.')
