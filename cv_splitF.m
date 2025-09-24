% Code adapted from: https://github.com/ThomasYeoLab/Ooi2022_MMP_HCP

%Define command line arguments.
function [testout] = cv_splitF(sublist,famlist,nfold,cseed)

%Get length.
nsub = size(sublist,2);

%Create family structure.
unique_famlist = unique(famlist); %Get unique family IDs.
nfam = length(unique_famlist); %Get number of unique families.
sub_perfam = cell(nfam,1); %Split subjects into families.
for fidx=1:nfam
    subidx = strcmp(famlist, unique_famlist{fidx})==1; %Find indices of current family.
    sub_perfam{fidx} = sublist(subidx)'; %Get subjects corresponding to the indices.
end

%Split families into folds with random seed.
fold_list = cell(nfold,1); %Generate a cell for each fold to collect.
subfold_amount = ceil(nsub/nfold); %Divide number of subjects by folds to see approximately how many will go into each fold.
rng(cseed); %Set seed.
rand_famidx = randperm(nfam); %Produce a randomized list of family indices.
cfold = 0; %Initialize the current fold to zero.
for fidx = 1:nfam %For each family.
    cfold = mod(cfold,nfold) + 1; %Find the remainder after dividing the current fold by the number of folds and add 1. Set this as the fold to put families in after incrementing.
    while size(fold_list{cfold},1)>=subfold_amount %If the number of subjects in this fold is more than the amount going into each fold, skip this fold and do the same incrementing.
        cfold = mod(cfold,nfold) + 1; 
    end
    fold_list{cfold} = [fold_list{cfold};sub_perfam{rand_famidx(fidx)}]; %Put the subjects of this family into the fold.
end

%Generate a table containing the subject and which test fold they are in.
testout = table(zeros(nsub,1),'RowNames',sublist);
for fidx=1:nfold
    cfold = fold_list{fidx};
    testout(cfold,1) = {fidx};
end
end
