Versions:

Packages.txt - Lists Python, MATLAB, and R versions and relevant packages with corresponding versions.

----------
Files:

atlasMvec2.csv - Upper triangular matrix vectorization indices.

colenetlabels.txt - Network labels retrieved from: https://github.com/ColeLab/ColeAnticevicNetPartition

_random_search.py - Replace corresponding version in the himalaya Python package.

_sklearn_api.py - Replace corresponding version in the himalaya Python package.

gcentroids.m - Finds the centroid of vectors according to a clustering.

cv_splitF.m - Generates split indices for cross-validation while respecting family structure.

cv_permF.m - Generates permutation indices for permutation testing while respecting family structure.

threshold_arbmeasure.m - Threshold edges in a square matrix by an arbitrary measure.

averMeas.m - Calculates average controllability on a square matrix.

moduMeas.m - Calculates modal controllability on a square matrix.

Gramian.m - Calculates the controllability Gramian on a square matrix.

----------
Analysis:

1. Extract ICA-FIX resting-state fMRI data and motion for the Project 1 dFC analysis.

    param_maker.py - Read HCP1200 subjects with resting-state fMRI scans based on the directory structure extracted from DataLad.

    fmri_get.txt - For each resting-state subject ID, extract the associated 3T REST1-2 ICA-FIX surface-based fMRI data files from the DataLad repository.

    voxdemean.py - For each resting-state scan for each subject ID, subtract mean BOLD timeseries voxelwise.

    meants_do.txt - For each resting-state subject ID, use the 360 cortical Glasser atlas to average demeaned voxelwise fMRI data.

    motion_get.txt - For each resting-state subject ID, extract the associated 3T REST1-2 fMRI relative FD at each data point from the DataLad repository.

2. Select subjects for the Project 1 dFC analysis.

    subselector.py - Read in all participants with ICA-FIX 3T resting-state data, then isolate participants who have all four scans. Use the HCP1200 data dictionary to isolate participants who have all 10 cognitive test scores of interest. Retrieve the mean relative FD for each participant for each run and exclude participants with higher than 0.2 in any run. Find the intersection across subjects.

3. Extract cognitive scores and confounds.

    cogscorer_PCA.py - Do PCA to extract the first PC score and loadings to use as gPCA. 

    subselector_FA.py - Relax the fMRI requirement to increase sample size to improve cognitive score extraction.  Get cognitive score data for expanded sample and add gPCA to it.

    cogscorer_FA.R - For the expanded sample, conduct EFA and CFA and retrieve scores and loadings.

    cogcon_prep.R - Gather cognitive scores and confounds into one file.

4. Do LEiDA.

    LE_dFC.m - For the mean ROI fMRI timeseries for each run and each participant, conduct LEiDA.

    LE_group_statecalc.m - For the LE(t) for each run and each participant, stack them and conduct k-medians clustering with a specified k for one random iteration of 500. 

    LE_group_iter_compare.py - Find the best iteration of k-medians clustering for each k. Values with lower within-cluster sum of distance (WCS) scores are better.

    LE_group_stateanalyze.m - For the best clustering, generate 2D matrices for LE(t) centroids and dFC(t) centroids by averaging according to the clustering. 

    LE_group_statesort.py - For the set of LE or dFC centroids across a range of k, print out the 2D images and generate a table which sorts them. Then, visually sort them and edit the table. Run this code again with the edited table to visualize the sorted plots.

    LE_group_stateclass.py - For the selected k, use the sorting matrix to translate the original clustering to the universal indices across k.

5. Generate dFC median states, SC, and sFC for each subject and thresholding inputs. 

    LE_group_substate.m - For each participant, find the median of the dFC(t) for each cluster for the clustering using the given k. 

    dFC_SC_sFC_analyze.m - Find the participants in Project 1 who also have SC data and generate a list. Then, for a k for dFC clustering and a type of SC normalization, generate the subject SC and sFC matrices. Then, find the coefficient of variation (CV) and median of these matrices to threshold the matrices by variability and strength later.

6. Calculate degree, average controllability, and modal controllability on the SC, dFC, and sFC for each subject.

    SC_subcontrol.m - For each participant, for a type of SC normalization, and the type and percentage of thresholding, generate the thresholded SC matrix, degree, average controllability, and modal controllability.

    SC_subcontrol_collect.m - Across participants, for a type of SC normalization, and the type and percentage of thresholding, collect the degree, average controllability, and modal controllability for SC.

    LE_group_subcontrol_structfunc.m - For each participant, for a k for clustering, and the type and percentage of thresholding, generate the thresholded dFC state matrices, degree, average controllability, and modal controllability.

    LE_group_subcontrol_structfunc_collect.m - For a k for clustering, and the type and percentage of thresholding, for each thresholded dFC state collect degree, average controllability, and modal controllability across participants.

    sFC_subcontrol_structfunc.m - For each participant, for the type and percentage of thresholding, generate the thresholded sFC matrix, degree, average controllability, and modal controllability.

    sFC_subcontrol_structfunc_collect.m -  Across participants, for the type and percentage of thresholding, collect the degree, average controllability, and modal controllability for sFC.

7. Calculate the principal cortical gradient from average sFC.

    sFC_subcontrol_structfunc_gradient.m - Across participants, read in their sFC matrices, Fisher Z-transform, average, reverse transform, and calculate gradients. Save the gradients, and check the variance and first three gradients.

    sFC_subcontrol_structfunc_gradient.R - For the given type and percentage of thresholding, plot each gradient for sFC on the brain. Generate a file to adjust the sign of the gradient as desired.

8. Calculate the controllability Gramians for the average dFC, sFC, and SC to remove regions with errors from further analyses. 
    
    LE_group_subcontrol_structfunc_stateview.m - Across participants, for a type of SC normalization, and the type and percentage of thresholding, average dFC, sFC, and SC matrices across participants and generate thresholded versions. Save the Gramian eigenvalues, average controllability, and modal controllability for each matrix. Also save which regions have errors regarding Gramian eigenvalues, to remove in later controllability analyses.

    LE_group_subcontrol_structfunc_stateview_controlqc.m - For a type of SC normalization, and the type and percentage of thresholding, for the average dFC, sFC, and SC matrices across participants, find the regions with errors regarding Gramian eigenvalues, verify the error in their corresponding values by collecting the Gramian eigenvalue, average controllability, and modal controllability. 

9. Plot the average matrices for SC, dFC, and sFC and the corresponding average and modal controllability on the brain.

    LE_group_subcontrol_structfunc_stateview.py - For the given k for clustering, SC type of normalization, and thresholding type and percentage, plot the dFC, SC, and sFC matrices organized according to hemisphere and network. Also save the network legends and color bars.

    LE_group_subcontrol_structfunc_stateview.R - For the given k for clustering, type of SC normalization, and type and percentage of thresholding, plot average controllability, modal controllability, and degree on the brain along with radar plots to examine regions with the highest values.

10. Find the correlations between cognitive variables, and between controllability and degree.

    LE_group_subcontrol_structfunc_allcorr.py - Find the correlations between g, gF, gC, gF-g, and gC-g for the main estimate with CFA. Also find the correlations for estimates based on EFA, PCA, and the NIH toolbox.

    LE_group_subcontrol_structfunc_allcorr.R - For the given k for clustering, type of SC normalization, and type and percentage of thresholding, generate box plots for the correlation between average controllability, modal controllability, and degree across subjects and across regions.

11. Generate Kernel Ridge Regression (KRR) models and comparisons of the models and do statistical tests. 

    LE_group_subcontrol_structfunc_predictF.m - For the given number of repetitions, outer loops, and inner loops for cross validation (CV) and the number of permutations, generate indices respecting family structure. Specifically, make sure families are not split between folds and make sure permutation only happens in blocks of unrelated participants.

    LE_group_subcontrol_structfunc_KRR_score.py - For a given k for clustering, type of SC normalization, and type and percentage of thresholding, fit a model with CV. This model predicts each cognitive variable in the given batch using either average controllability, modal controllability, or degree calculated from the set of state matrices specified. CV is done for the given repetitions, and inner and outer fold number. To facilitate permutation testing with the same parameters set, for each fold, the Y variable with transformations applied prior to model fitting, alpha, beta, and gamma are all saved. Next, the CV accuracy for the test set, the CV accuracy for the inner CV for the best hyperparameter set, the decomposed CV accuracy for the test set for each feature band, and the decomposed CV accuracy after demeaning of the test set according to the decomposition assumptions are all saved. Lastly, the Haufe scores are also saved.

    LE_group_subcontrol_structfunc_KRR_score_format.R - For a given k for clustering, type of SC normalization, and type and percentage of thresholding, collect decomposed accuracy scores on the model predicting each cognitive variable in the given batch using either average controllability, modal controllability, or degree calculated from the set of state matrices specified. Reformat to generate percentages of the total accuracy score.

    LE_group_subcontrol_structfunc_KRR_perm.py - For a given k for clustering, type of SC normalization, and type and percentage of thresholding, fit models from CV with a single permutation using the same parameters from the original model. The permuted CV accuracy for the test set and Haufe scores are both saved.

    LE_group_subcontrol_structfunc_KRR_analyze.py - For a given k for clustering, type of SC normalization, and type and percentage of thresholding, do permutation testing on the model predicting each cognitive variable in the given batch using either average controllability, modal controllability, or degree calculated from the set of state matrices specified. The one-sided p-value is produced from the test set CV accuracy score because the Coefficient of Determination (COD) is only positive, and two-sided p-values are produced for the Haufe scores because Haufe scores can be positive or negative. Accuracy, decomposed accuracy, and Haufe scores are all packaged for retrieval.

    LE_group_subcontrol_structfunc_KRR_compare.py - For a given k for clustering, type of SC normalization, and type and percentage of thresholding, do corrected paired t-tests between the different models predicting each cognitive variable using either average controllability, modal controllability, or degree calculated from the set of state matrices specified.

    LE_group_subcontrol_structfunc_KRR_scorcomp_FDR.m - Across participants, for a type of SC normalization, and the type and percentage of thresholding, do FDR across the test set CV accuracy scores and accuracy comparison scores of interest.

    LE_group_subcontrol_structfunc_KRR_feat_FDR.m - Across participants, for a type of SC normalization, and the type and percentage of thresholding, do FDR across the Haufe feature scores from all models of interest.

    LE_group_subcontrol_structfunc_KRR_detailed_plot.R - For a given k for clustering, type of SC normalization, type and percentage of thresholding, controllability type, the set of matrices, the set of cognitive variables, CV repetitions, CV inner fold number, CV outer fold number, and permutation number, plot the regions with significant Haufe scores and corresponding radar plots for the networks and regions to analyze.

    LE_group_subcontrol_structfunc_KRR_common_plot.R - For a given k for clustering, type of SC normalization, type and percentage of thresholding, the set of controllability, the set of matrices, the set of cognitive variables, CV repetitions, CV inner fold number, CV outer fold number, and permutation number, plot the regions with significant Haufe scores which are common and unique between average and modal controllability.

12. Investigate correlations between Haufe scores and different brain maps.

    LE_group_subcontrol_structfunc_KRR_maplook.py - For a given k for clustering, type of SC normalization, and type and percentage of thresholding, find Dice coefficients between significant Haufe scores and PFIT/MD maps. Generate p-values for Dice coefficients using permutation. Find Spearman's correlation to quantify the relationship and R2 for quadratic regression to quantify the statistical relationship between Haufe scores and controllability/principal cortical gradient. Generate p-values for the R2.

    LE_group_subcontrol_structfunc_KRR_maplook_plot.py - For a given k for clustering, type of SC normalization, and type and percentage of thresholding, illustrate the relationship between Haufe scores and controllability/principal cortical gradient using linear and quadratic regression lines and Pearson's and Spearman's correlation.

    LE_group_subcontrol_structfunc_KRR_maplook_FDR.m - Across participants, for a type of SC normalization, and the type and percentage of thresholding, retrieve the brain map analysis statistically tested values of interest to FDR correct. These include the Dice coefficients between significant Haufe scores and PFIT/MD maps and R2 for quadratic regression to quantify the statistical relationship between Haufe scores and controllability/principal cortical gradient. 

13. Plot the accuracy scores for models while varying parameters.

    LE_group_subcontrol_structfunc_KRR_paramvar.py - For a given number of repetitions, inner loop number, and outer loop number for CV, collect the CV average test set accuracy score from the models using different variations of average or modal controllability, matrix modalities, thresholding types, thresholding values, and cognitive measures. First build the labels, then using the labels, read them all into a matrix for plotting use.

    LE_group_subcontrol_structfunc_KRR_paramvar_plot.R - For a given number of repetitions, inner loop number, and outer loop number for CV, plot the CV average test set accuracy score for different models to illustrate changes from varying the threshold type and value, k for clustering, and how the cognitive variable is measured.

14. Investigate what happens when both average and modal controllability are input in the same model.

    LE_group_subcontrol_structfunc_KRR_score_bothctrl.py - For a given k for clustering, type of SC normalization, and type and percentage of thresholding, fit a model with CV similar to the main analysis but with both average and modal controllability in the same model. Generate the CV accuracy for the test set, the CV accuracy for the inner CV for the best hyperparameter set, the decomposed CV accuracy for the test set for each feature band, the decomposed CV accuracy after demeaning of the test set according to the decomposition assumptions, and the Haufe scores.

    LE_group_subcontrol_structfunc_KRR_score_format_bothctrl.R - For a given k for clustering, type of SC normalization, and type and percentage of thresholding, collect decomposed accuracy scores on the model predicting each cognitive variable in the given batch using both average and modal controllability calculated from the set of state matrices specified. Reformat to generate percentages of the total accuracy score.
