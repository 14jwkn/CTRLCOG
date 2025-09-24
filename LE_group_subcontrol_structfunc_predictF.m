%{
For the given number of repetitions, outer loops, and inner loops for cross
validation (CV) and the number of permutations, generate indices respecting
family structure. Specifically, make sure families are not split between
folds and make sure permutation only happens in blocks of unrelated
participants.
Outputs:
'r',num2str(nrep),'_o',num2str(outer_k),'_i',num2str(inner_k),'_predict_splits.h5'
Outer and inner cross validation subject indices.
'_perm_splits.h5' Permutation subject indices.
%}

%Define command line arguments.
function [] = LE_group_subcontrol_structfunc_predictF(nrep,inner_k,outer_k,nperm)
% nrep = Number of repetitions for the CV.
% inner_k = Number of inner loops for the CV.
% outer_k = Number of outer loops for the CV.
% nperm = Number of permutations for permutation testing.
disp(append('Doing: ',nrep,' ',inner_k,' ',outer_k,' ',nperm))
        
%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))        
        
%Generate output path and filename.
subgroup = 'full';
sc_subgroup = 'dr_full';
outpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                  subgroup,'/',k,'/SC_dFC/',sc_subgroup,'/predict_splits/');
              
%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end            
        
%Convert to numeric.
nrep = str2num(nrep);
inner_k = str2num(inner_k);
outer_k = str2num(outer_k);
nperm = str2num(nperm);
        
%Read in the subjects and family structure.
infile = 'dr_full_intersect.txt';
sublist = cellstr(string(readmatrix(infile)))';
nsub = size(sublist,2);
infile = '../inputs/data/hcp/RESTRICTED_HCP1200_Data_Dictionary.csv';
demomat = readtable(infile,ReadRowNames=true);
famlist = table2cell(demomat(sublist,'Family_ID'))';

%Set output file for CV.
outfile = append(outpath,'r',num2str(nrep),'_o',num2str(outer_k),'_i',num2str(inner_k),'_predict_splits.h5');

%Initialize starting seed for rng.
seedstart = 12345;

%Collect repetitions.
rep_testcollect = zeros(nsub,nrep);

%For each repetition.
repcount_want = 1;
repcount_true = 1;
while repcount_want < (nrep+1)
   
    %Initialize seed and do splitting for outer CV.
    disp(repcount_true)
    cseed = repcount_true + seedstart - 1;
    testout = cv_splitF(sublist,famlist,outer_k,cseed);
    
    %If this split is not done already.
    testarr = table2array(testout);
    testcomp = ismember(testarr',rep_testcollect','rows');
    if ~testcomp
 
        %Put the splits into the matrix.
        rep_testcollect(:,repcount_want) = testarr;
        
        %Increment number of repetitions wanted.
        repcount_want = repcount_want + 1;
    end
    
    %Increment true number of repetitions.
    repcount_true = repcount_true + 1;
end

%Save the outer CV folds for each repetition.
outkey = append('/outercv_testidx');
try
    h5create(outfile,outkey,size(rep_testcollect));
catch
    disp('File and key created already.')
end
h5write(outfile,outkey,rep_testcollect);

%For each repetition.
for ridx=1:nrep
    
    %Extract outer fold test indices for this repetition.
    cout = rep_testcollect(:,ridx);
    
    %Collect outer folds.
    outer_testcollect = zeros(nsub,outer_k);
 
    %For each outer fold.
    for ouidx=1:outer_k
        
        %Extract outer fold train indices for this outer fold.
        cout_train = (cout ~= ouidx);
        
        %Generate subject list and family list for these outer fold train indices.
        sublist_train = sublist(cout_train);
        famlist_train = famlist(cout_train);

        %Initialize seed and do splitting for inner CV.
        cseed = ouidx + seedstart - 1;
        testout = cv_splitF(sublist_train,famlist_train,inner_k,cseed);
        
        %Put the splits into the matrix.
        outer_testcollect(cout_train,ouidx) = table2array(testout);
    end
    
    %Save the inner CV folds for each outer CV fold.
    outkey = append('/r',num2str(ridx),'_innercv_testidx');
    try
        h5create(outfile,outkey,size(outer_testcollect));
    catch
        disp('File and key created already.')
    end
    h5write(outfile,outkey,outer_testcollect); 
end

%Set output file for permutation.
outfile = append(outpath,num2str(nperm),'_perm_splits.h5');

%Collect permutations.
permcollect = zeros(nsub,nperm);

%For each permutation.
permcount_want = 1;
permcount_true = 1;
while permcount_want < (nperm+1)
   
    %Initialize seed and do permutation.
    cseed = permcount_true + seedstart - 1;
    permarr = cv_permF(sublist,famlist,cseed);
    
    %If this permutation is not done already.
    permcomp = ismember(permarr',permcollect','rows');
    if ~permcomp
 
        %Put the permutation into the matrix.
        permcollect(:,permcount_want) = permarr;
        
        %Increment number of permutations wanted.
        permcount_want = permcount_want + 1;
    end
    
    %Increment true number of permutations.
    permcount_true = permcount_true + 1;
end

%Save the permutations.
outkey = append('/permidx');
try
    h5create(outfile,outkey,size(permcollect));
catch
    disp('File and key created already.')
end
h5write(outfile,outkey,permcollect);
disp('Saved.')
end
    