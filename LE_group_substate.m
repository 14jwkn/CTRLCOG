%{
For each participant, find the median of the dFC(t) for each cluster for the
clustering using the given k. 
Output:
med_dFC.h5 Median of the dFC states for the participant.
%}

%Define command line arguments.
function [] = LE_group_substate(k,subject)
% k = Selected k for k-medians clustering.
% subject = Subject ID.
disp(append('Doing: ',k,' ',subject))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O. 
subfile = 'r_full_submain.txt';
subgroup = 'full'; 
outpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                 subgroup,'/',k,'/',subject,'/');
dFCmed_file = append(outpath,'med_dFC.h5');

%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end

%Read in subjects.
subjects = textread(subfile,'%s','delimiter',',');
subjects = string(subjects); 

%Get the subject index and produce the start and end indices.
subin = find(strcmp(subjects,subject));
if ~strcmp(order,'both')
    nwin = 2396;
else
    nwin = 4792;
end
startin = (nwin*(subin-1)) + 1;
endin = nwin*subin;
disp(append('Index: ',num2str(startin),'-',num2str(endin)))

%Read in the best iteration.
best_file = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/',...
                  'group/best_iter/',subgroup,'/best_iter.csv');
best_clust = readtable(best_file);
best_lab = k;
iteration = num2str(best_clust{:,best_lab});
disp(append('Best: ',iteration))

%Read in best clustering.
inpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                subgroup,'/',k,'/');
clustfile = append(inpath,'subclust_',iteration,'.h5');
inkey = append('/',subgroup);
bestclust = h5read(clustfile,inkey);
bestclust = bestclust(startin:endin);

%Read in the dFC data.
Lfile1 = append('../outputs/r_LE_dFC/',clean,'/REST1_LR/',subject,...
                '/LE_dFC.h5');     
Rfile1 = append('../outputs/r_LE_dFC/',clean,'/REST1_RL/',subject,...
                '/LE_dFC.h5'); 
Lfile2 = append('../outputs/r_LE_dFC/',clean,'/REST2_LR/',subject,...
                '/LE_dFC.h5'); 
Rfile2 = append('../outputs/r_LE_dFC/',clean,'/REST2_RL/',subject,...
                 '/LE_dFC.h5');  

%Index the h5.
inkey = '/LE_dFC';
dFCmatL1 = h5read(Lfile1,inkey);
dFCmatR1 = h5read(Rfile1,inkey);
dFCmatL2 = h5read(Lfile2,inkey);
dFCmatR2 = h5read(Rfile2,inkey);
dFCmat = vertcat(dFCmatL1,dFCmatR1,dFCmatL2,dFCmatR2);
delvars = {'dFCmatL1','dFCmatR1','dFCmatL2','dFCmatR2'};
clear(delvars{:})  
disp('Read input matrices.')

%Produce the mean centroids.
clusts = 1:str2num(k);
distance = 'cityblock';
[dFCmed,~] = gcentroids(dFCmat,bestclust,clusts,distance);      

%Save as tables.
outkey = append('/',subgroup);
try
    h5create(dFCmed_file,outkey,size(dFCmed));
catch
    disp('File and key created already.')
end
h5write(dFCmed_file,outkey,dFCmed);
disp('Saved.')
end
