% Code adapted from: https://github.com/ThomasYeoLab/Ooi2022_MMP_HCP

%Define command line arguments.
function [y_perm] = cv_permF(sublist,famlist,cseed)

%Set up labels.
subvec = string(sublist); 
nsub = size(subvec,2);
famvec = string(famlist);

%Extract convert family IDs to family indices, get number of unique families,
%and get the size of each family.
[~, ~,famidx] = unique(famvec); 
unique_nfam = max(famidx);
famsize = zeros(unique_nfam,1);
for fidx = 1:unique_nfam
	famsize(fidx,:) = sum(famidx==fidx);
end

%Initialize seed.
rng(cseed)

%Collect permutation groups containing different families for each subject.
perm_group = zeros(nsub,1);

%For each family.
for ufidx = 1:unique_nfam
    
    %For the family, get the count of subjects. From 1 to the max number of
    %subjects per family, get a number for each subject in the family.
    %Basically randomly splitting the subjects of a family into different groups,
    %so that permutations do not occur in groups containing the same
    %family.
    sorting_order = datasample(1:max(famsize),famsize(ufidx),'Replace',false);
    
    %Record the number for the subject.
    perm_group(famidx==ufidx) = sorting_order; 
end

%Set up indices and permuted indices.
y_idx = (1:nsub)';
y_perm = zeros(nsub,1); 

%For each group, permute indices within the groups.
for groupidx = 1:max(famsize)
    y_tmp = y_idx(perm_group == groupidx);
    y_tmp = y_tmp(randperm(length(y_tmp)));
    y_perm(perm_group == groupidx) = y_tmp;
end
end
