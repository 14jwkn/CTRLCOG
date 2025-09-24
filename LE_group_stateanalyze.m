%{
For the best clustering, generate 2D matrices for LE(t) and dFC(t) centroids by 
averaging according to  the clustering. Next, average while multiplying the 
dFC(t) by occurrence.
Adapts code from: 
https://github.com/juanitacabral/LEiDA
https://github.com/trendscenter/gift
Output:
k_LE.h5 2D centroids for LE(t).
k_dFC.h5 2D centroids for dFC(t).
%}

%Define command line arguments.
function [] = LE_group_stateanalyze(k)
disp(append('Doing ',k));

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O. 
subfile = 'r_full_submain.txt';
subgroup = 'full';       
outpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                subgroup,'/',k,'/');

%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end

%Read in subjects.
subjects = textread(subfile,'%s','delimiter',',');
subjects = string(subjects); 
nsubj = size(subjects,1);

%Read in the best iteration.
best_file = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/',...
                  'group/best_iter/',subgroup,'/best_iter.csv');
best_clust = readtable(best_file);
best_lab = k;
iteration = num2str(best_clust{:,best_lab});
disp(append('Best: ',iteration))

%Read in best clustering and centroids, then initialize values.
clustfile = append(outpath,'subclust_',iteration,'.h5');
centfile = append(outpath,'subcent_',iteration,'.h5');
inkey = append('/',subgroup);
IDX = h5read(clustfile,inkey);
V = h5read(centfile,inkey);
Tmax = size(IDX,1);
[N_Cl,N_ba] = size(V);

%Do LE.
disp('Doing LE')
tic
outfile = append(outpath,'k_LE.h5');
for c=1:N_Cl
    
    %Produce 2D matrix.
    VVT=V(c,:)'*V(c,:); 
    
    %Save 2D matrix.
    outkey = append('/',subgroup,'_',num2str(c));
    try
        h5create(outfile,outkey,size(VVT));
    catch
        disp('File and key created already.')
    end
    h5write(outfile,outkey,VVT);
end
disp('LE saved.')
toc

%Do dFC.
tic
nconn = (N_ba*(N_ba-1))/2;
k_mean = zeros(N_Cl,nconn);
for subin = 1:nsubj
        
    %Read in subject data.
    disp(append('Doing dFC: ',subjects(subin)))
    inkey = '/LE_dFC';
    Lfile1 = append('../outputs/r_LE_dFC/REST1_LR/',subjects(subin),...
                        '/LE_dFC.h5');     
    Rfile1 = append('../outputs/r_LE_dFC/REST1_RL/',subjects(subin),...
                    '/LE_dFC.h5'); 
    Lfile2 = append('../outputs/r_LE_dFC/REST2_LR/',subjects(subin),...
                        '/LE_dFC.h5'); 
    Rfile2 = append('../outputs/r_LE_dFC/REST2_RL/',subjects(subin),...
                    '/LE_dFC.h5'); 
    dFCmatL1 = h5read(Lfile1,inkey);
    dFCmatR1 = h5read(Rfile1,inkey);
    dFCmatL2 = h5read(Lfile2,inkey);
    dFCmatR2 = h5read(Rfile2,inkey);
    dFCmat = vertcat(dFCmatL1,dFCmatR1,dFCmatL2,dFCmatR2);
    delvars = {'dFCmatL1','dFCmatR1','dFCmatL2','dFCmatR2'};
    clear(delvars{:})  
    nwin = size(dFCmat,1);

    %Set time indices and go through them.
    startin = (nwin*(subin-1)) + 1;
    endin = nwin*subin;
    subt = 1;
    disp(append('Subject t: ',num2str(startin),' to ',num2str(endin)));
    for t=startin:endin

        %Add to the corresponding running k unweighted mean.
        idx_c = IDX(t);
        k_mean(idx_c,:) = k_mean(idx_c,:) + (dFCmat(subt,:));

        %Increment subject time.
        subt = subt + 1; 
    end
    delvars = {'dFCmat'};
    clear(delvars{:})
end
k_mean = k_mean/Tmax;

%Unvectorize and save the k means.
outfile = append(outpath,'k_dFC.h5');
for c=1:N_Cl
    
    %Produce 2D matrix. 
    k_2D = zeros(N_ba,N_ba);
    k_2D(triu(true(N_ba),1)) = k_mean(c,:);
    k_2D = k_2D + k_2D.';
    
    %Save 2D matrix.
    outkey = append('/',subgroup,'_',num2str(c));
    try
        h5create(outfile,outkey,size(k_2D));
    catch
        disp('File and key created already.')
    end
    h5write(outfile,outkey,k_2D);
end
disp('dFC saved.')
toc
end
