%{
For the LE(t) for each run and each participant, stack them and conduct k-medians 
clustering with a specified k for one random iteration of 500. 
Output:
subclust_*.h5 Clustering for this iteration.
subcent_*.h5 h5 Centroids for this iteration.
substat_*.h5 Clustering-based metrics for this iteration.
%}

%Define command line arguments.
function [] = LE_group_statecalc(k,iteration)
%k = Selected k for k-medians clustering.
%iteration = Current random iteration of 500.
disp(append('Doing ',k,' ',iteration));
tic

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O.
subfile = 'r_full_submain.txt';
subgroup = 'full';      
outpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
    subgroup,'/',k,'/');
clustfile = append(outpath,'subclust_',iteration,'.h5');
centfile = append(outpath,'subcent_',iteration,'.h5');
statfile = append(outpath,'substat_',iteration,'.h5');

%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end 

%Read in subjects.
subjects = textread(subfile,'%s','delimiter',',');
subjects = string(subjects); 
nsubj = size(subjects,1);

%Read in first subject.
inkey = '/LE_LE';
Lfile1 = append('../outputs/r_LE_dFC/REST1_LR/',subjects(1),...
                    '/LE_dFC.h5');     
Rfile1 = append('../outputs/r_LE_dFC/REST1_RL/',subjects(1),...
                '/LE_dFC.h5'); 
Lfile2 = append('../outputs/r_LE_dFC/REST2_LR/',subjects(1),...
                    '/LE_dFC.h5'); 
Rfile2 = append('../outputs/r_LE_dFC/REST2_RL/',subjects(1),...
                '/LE_dFC.h5'); 
dFCmatL1 = h5read(Lfile1,inkey);
dFCmatL1 = dFCmatL1(:,roi_index);
dFCmatR1 = h5read(Rfile1,inkey);
dFCmatR1 = dFCmatR1(:,roi_index);
dFCmatL2 = h5read(Lfile2,inkey);
dFCmatL2 = dFCmatL2(:,roi_index);
dFCmatR2 = h5read(Rfile2,inkey);
dFCmatR2 = dFCmatR2(:,roi_index);
dFCmat = vertcat(dFCmatL1,dFCmatR1,dFCmatL2,dFCmatR2);
delvars = {'dFCmatL1','dFCmatR1','dFCmatL2','dFCmatR2'};
clear(delvars{:})  
nwin = size(dFCmat,1);
nconn = size(dFCmat,2);

%Set up empty matrix and insert first subject's data.
dFCfull = zeros(nwin*nsubj,nconn);
dFCfull(1:nwin,:) = dFCmat;
delvars = {'dFCmat'};
clear(delvars{:})

%For each subject.
for subin = 2:nsubj
    
    %Set indices.
    startin = (nwin*(subin-1)) + 1;
    endin = nwin*subin;
    
    %Read in and append. 
    Lfile1 = append('../outputs/r_LE_dFC/REST1_LR/',subjects(subin),...
                        '/LE_dFC.h5');     
    Rfile1 = append('../outputs/r_LE_dFC/REST1_RL/',subjects(subin),...
                    '/LE_dFC.h5'); 
    Lfile2 = append('../outputs/r_LE_dFC/REST2_LR/',subjects(subin),...
                        '/LE_dFC.h5'); 
    Rfile2 = append('../outputs/r_LE_dFC/REST2_RL/',subjects(subin),...
                    '/LE_dFC.h5'); 
    dFCmatL1 = h5read(Lfile1,inkey);
    dFCmatL1 = dFCmatL1(:,roi_index);
    dFCmatR1 = h5read(Rfile1,inkey);
    dFCmatR1 = dFCmatR1(:,roi_index);
    dFCmatL2 = h5read(Lfile2,inkey);
    dFCmatL2 = dFCmatL2(:,roi_index);
    dFCmatR2 = h5read(Rfile2,inkey);
    dFCmatR2 = dFCmatR2(:,roi_index);
    dFCfull(startin:endin,:) = vertcat(dFCmatL1,dFCmatR1,dFCmatL2,dFCmatR2);
    delvars = {'dFCmatL1','dFCmatR1','dFCmatL2','dFCmatR2'};
    clear(delvars{:})   
end
disp('Read dFC matrix.')
toc

%Set up parameters.
num_clusters = str2num(k);
dmethod = 'cityblock';
kmeans_num_replicates = 1; 
kmeans_max_iter = 10000;
dispmethod = 'final';
emptytreat = 'drop';
randseed = 12345*str2num(iteration);

%Run k-clustering.
tic
rng(randseed)
[IDXp, Cp, SUMDp, ~] = kmeans(dFCfull,num_clusters,'distance', dmethod,...
    'Replicates',kmeans_num_replicates,'MaxIter',kmeans_max_iter,...
    'Display',dispmethod,'empty',emptytreat);
disp('Finished clustering.')
toc

%Package the output.
packout = SUMDp.';

%Try making the files, if it already exists then move on.
outkey = append('/',subgroup);
try
    h5create(clustfile,outkey,size(IDXp));
catch
    disp('File and key created already.')
end
try
    h5create(centfile,outkey,size(Cp));
catch
    disp('File and key created already.')
end
try
    h5create(statfile,outkey,size(packout));
catch
    disp('File and key created already.')
end

%Write out outputs.
h5write(clustfile,outkey,IDXp);
h5write(centfile,outkey,Cp);
h5write(statfile,outkey,packout);
disp('Saved.')
end
