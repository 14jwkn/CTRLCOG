%{
For the mean ROI fMRI timeseries for each run and each participant, conduct 
LEiDA. Adapts code from: 
https://github.com/juanitacabral/LEiDA
Output:
LE_dFC.h5 Phase dFC(t) and LE(t) for all time points. 
%}

%Define command line arguments.
function [] = LE_dFC(run,subject)
%run = run_label
%subject = subject_ID
disp(append('Doing ',run,' ',subject));

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O.
inpath = append('../outputs/r_meants/',subject,'/');
infile = append('demean_rfMRI_',run,'_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv');
outpath = append('../outputs/r_LE_dFC/',run,'/',subject,'/')  ; 
outfile = append(outpath,'LE_dFC.h5'); 

%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end

%Set up parameters.
tc = readmatrix(append(inpath,infile))';
[all_ts, N_areas] = size(tc);
Tmax = all_ts - 2;

%Set up empty dFC, LE, and LE variance matrices.
dFC=zeros(Tmax,N_areas,N_areas); 
Leading_Eig=zeros(Tmax,N_areas); 
    
%Get the BOLD phase using the Hilbert transform.
BOLD = tc.';
Phase_BOLD=zeros(N_areas,all_ts);
for seed=1:N_areas
    BOLD(seed,:) = BOLD(seed,:)-mean(BOLD(seed,:)); %Demean.
    Phase_BOLD(seed,:) = angle(hilbert(BOLD(seed,:))); %Hilbert transform.
end
clear tc

%Get rid of the last and first iFC TRs for signal distortions.
Phase_BOLD(:,all_ts) = [];
Phase_BOLD(:,1) = [];

%For each time point.
for t=1:Tmax
    
    %Calculate the instantaneous phase coherence, iFC.
    iFC=zeros(N_areas);
    for n=1:N_areas
        for p=1:N_areas
            iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
        end
    end
    dFC(t,:,:) = iFC;
    
    %Calculate the leading eigenvector.
    [V1,~]=eigs(iFC,1);
    if mean(V1>0)>.5
        V1=-V1;
    elseif mean(V1>0)==.5 && sum(V1(V1>0))>-sum(V1(V1<0))
        V1=-V1;
    end
    Leading_Eig(t,:)=V1;
end 

%Vectorize dFC matrix.
I_sup_diag=find(triu(ones(N_areas),1));
nconn=size(I_sup_diag,1);
dFC_vec = zeros(Tmax,nconn);
for t=1:Tmax
    dFC_mat = squeeze(dFC(t,:,:));
    dFC_vec(t,:) = dFC_mat(I_sup_diag);
end    

%Try making the file, if it already exists then move on.
dFCkey = '/LE_dFC';
lekey = '/LE_LE';
try
    h5create(outfile,dFCkey,size(dFC_vec));
catch
    disp('File and key created already.')
end
try
    h5create(outfile,lekey,size(Leading_Eig));
catch
    disp('File and key created already.')
end

%Write out dFC files.
h5write(outfile,dFCkey,dFC_vec);
h5write(outfile,lekey,Leading_Eig);
disp('Saved.')
end
