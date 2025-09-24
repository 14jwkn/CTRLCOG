%{
Find the participants in Project 1 who also have SC data and generate a list. 
Then, for a k for dFC clustering and a type of SC normalization, generate the 
subject SC and sFC matrices. Then, find the coefficient of variation (CV) 
and median of these matrices to threshold the matrices by variability and strength later.
Output:
dr_full_intersect.txt List of participants in Project 1 with SC data.
'sub_sc_',sctype,'.h5' SC matrix for each participant.
sub_sFC.h5 sFC matrix for each participant.
'group_cv_',sctype,'.csv' Group CV for SC.
'group_med_',sctype,'.csv' Group median for SC.
group_cv_dFC.csv Group CV for dFC.
group_med_dFC.csv Group median for dFC.
group_cv_sFC.csv Group CV for sFC.
group_med_sFC.csv Group median for sFC.
%}

%Define command line arguments.
function [] = dFC_SC_sFC_analyze(k,sctype)
% k = Selected k for k-medians clustering.
% sctype = Type of SC normalization.
disp(append('Doing: ',k,' ',sctype))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set parameters.
nk = str2num(k);
nroi = 360;
nconn = 64620;
thresvals = [0 10 20 30 40 50 60 70 80 90];
nthres = size(thresvals,2);
vecmaker = find(triu(ones(nroi),1));

%Set outpath.
subgroup = 'full';
sc_subgroup = 'dr_full';
outpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/', ...
                 subgroup,'/',k,'/SC_dFC/',sc_subgroup,'/');
if not(isfolder(outpath))
    mkdir(outpath)
end

%Read in subjects selecting for resting-state FC data.
subfile = 'r_full_submain.txt';       
r_subjects = textread(subfile,'%s','delimiter',',');
r_subjects = string(r_subjects); 

%Get centroid resorting.
classfile = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                       'visclass/',subgroup,'/allk/',subgroup,'_sortclass.csv');
classmat = readtable(classfile,'ReadRowNames',true,'ReadVariableNames',true);
classmat.Properties.VariableNames = string(1:size(classmat,2));
classmat = classmat(k,:);
classmat = rmmissing(classmat,2);

%Read in SC matrices. Download to the given path from https://zenodo.org/records/4060485.
if strcmp(sctype,'fpt') 
    infile = '../inputs/data/dti/individualConnectivity_10^Fpt.mat';
    groupd = load(infile);
elseif strcmp(sctype,'raw')
    infile = '../inputs/data/dti/individualConnectivity_rawStreamlineCount.mat';
    groupd = load(infile);
end 

%Extract subjects with SC and find matches with resting-state FC analysis subjects.
d_subjects = string(groupd.subjectIDs);
dr_subjects = intersect(r_subjects,d_subjects);
outfile = 'dr_full_intersect.txt';
writematrix(dr_subjects,outfile)

%Read in subject SC matrices.
wantsubs = dr_subjects;
nsub = size(wantsubs,1);
sc_stack = zeros(nsub,nconn);
for sidx=1:nsub
    
    %Extract.
    csub = wantsubs(sidx);
    cidx = find(d_subjects==csub);
    cmat = groupd.rawStreamlineCounts(:,:,cidx);
    
    %Resort parcels for our atlas order.
    our_order = [181:360,1:180];
    cmat = cmat(our_order,our_order);
    
    %Vectorize and add.
    sc_stack(sidx,:) = cmat(vecmaker);
    
    %Save.
    outkey = append('/',csub);  
    subpath = append('../inputs/data/dti/',csub,'/');
    if not(isfolder(subpath))
        mkdir(subpath)
    end
    outfile = append(subpath,'sub_sc_',sctype,'.h5');
    try
        h5create(outfile,outkey,size(cmat));
    catch
        disp('File and key created already.')
    end
    h5write(outfile,outkey,cmat);
end    

%Read in subject dFC state matrices.
state_stack = cell(1,nk);
for kidx=1:nk
    
    %Extract.
    ck = num2str(kidx);
    cstack = zeros(nsub,nconn);
    disp(ck)
    
    %Go through each subject.
    for sidx=1:nsub
        
        %Extract.
        csub = wantsubs(sidx);
        
        %Read in.
        inpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                        subgroup,'/',k,'/',csub,'/');
        infile = append(inpath,'med_dFC.h5'); 
        inkey = append('/',subgroup);
        subcent = h5read(infile,inkey);
        
        %Resort states.
        subcent = subcent(table2array(classmat),:);

        %Extract and append.
        cstack(sidx,:) = subcent(kidx,:);        
    end
    
    %Add to matrix.
    state_stack{kidx} = cstack;
end

%Produce subject sFC matrices.
sFC_stack = zeros(nsub,nconn);
sFC_stack_fish = zeros(nsub,nconn);
for sidx=1:nsub

    %Extract.
    csub = wantsubs(sidx);

    %Read in runs.
    inpath = append('../outputs/r_meants/',csub,'/');
    prefix = append('demean_rfMRI_');

    %Set files.
    Lfile1 = append(inpath,prefix,'REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv');
    Rfile1 = append(inpath,prefix,'REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv');              
    Lfile2 = append(inpath,prefix,'REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv');
    Rfile2 = append(inpath,prefix,'REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii_meants.csv');  
    
    %Read and remove first and last time points to match with LEiDA dFC.
    sub_tc_L1 = readmatrix(Lfile1)';
    sub_tc_L1(1200,:) = [];
	sub_tc_L1(1,:) = [];
    sub_tc_R1 = readmatrix(Rfile1)';
    sub_tc_R1(1200,:) = [];
	sub_tc_R1(1,:) = [];
    sub_tc_L2 = readmatrix(Lfile2)';
    sub_tc_L2(1200,:) = [];
	sub_tc_L2(1,:) = [];
    sub_tc_R2 = readmatrix(Rfile2)';
    sub_tc_R2(1200,:) = [];
	sub_tc_R2(1,:) = [];
    
    %Stack them, correlate, and r-to-z transform.
    subcat = vertcat(sub_tc_L1,sub_tc_R1,sub_tc_L2,sub_tc_R2);
    subcorr = corrcoef(subcat);
    subcorr = subcorr(vecmaker);
    subcorr_fish = atanh(subcorr);
    
    %Append.    
    sFC_stack(sidx,:) = subcorr;  
    sFC_stack_fish(sidx,:) = subcorr_fish;   
end

%Save the sFC stack.
sFCpath = append('../outputs/r_sFC/',sc_subgroup,'/');
if not(isfolder(sFCpath))
    mkdir(sFCpath)
end 
outfile = append(sFCpath,'sub_sFC.h5');
outkey = append('/sFC');
try
    h5create(outfile,outkey,size(sFC_stack));
catch
    disp('File and key created already.')
end
h5write(outfile,outkey,sFC_stack);

%Generate CV for each matrix.
intab = sc_stack;
avgtab = mean(intab);
stdtab = std(intab);
sc_cv = stdtab./avgtab;
state_cv = zeros(nk,nconn);
for kidx=1:nk
    cstack = state_stack{kidx};
    intab = cstack;
    avgtab = mean(intab);
    stdtab = std(intab);
    state_cv(kidx,:) = stdtab./avgtab;
end
intab = sFC_stack;
avgtab = mean(intab);
stdtab = std(intab);
sFC_cv = stdtab./avgtab;

%Generate SC medians. 
sc_med = median(sc_stack,1);

%Generate dFC state medians.
state_med = zeros(nk,nconn);
for kidx=1:nk
    
    %Extract.
    cstack = state_stack{kidx};
    ccent = median(cstack,1);
    state_med(kidx,:) = ccent;
end

%Generate sFC medians using the Fisher transformed values to reduce bias and
%convert back after.
sFC_med_fish = median(sFC_stack_fish,1);
sFC_med = tanh(sFC_med_fish);

%Save group CV for SC.
dtipath = append('../outputs/d_SC/',sc_subgroup,'/');
if not(isfolder(dtipath))
    mkdir(dtipath)
end
outmat = sc_cv;
outtab = array2table(outmat);     
outfile = append(dtipath,'group_cv_',sctype,'.csv');
writetable(outtab,outfile,'WriteRowNames',0,'WriteVariableNames',0) 

%Save group strength for SC.
outmat = sc_med;
outtab = array2table(outmat);     
outfile = append(dtipath,'group_med_',sctype,'.csv');
writetable(outtab,outfile,'WriteRowNames',0,'WriteVariableNames',0) 

%Save group CV for dFC states.
if strcmp(k,'11')
    outmat = state_cv';
else
    outmat = state_cv;
end
outtab = array2table(outmat);     
outfile = append(outpath,'group_cv_dFC.csv');
writetable(outtab,outfile,'WriteRowNames',0,'WriteVariableNames',0) 

%Save group strength for dFC states.
if strcmp(k,'11')
    outmat = state_med';
else
    outmat = state_med;
end
outtab = array2table(outmat);     
outfile = append(outpath,'group_med_dFC.csv');
writetable(outtab,outfile,'WriteRowNames',0,'WriteVariableNames',0) 

%Save group CV for sFC.
sFCpath = append('../outputs/r_sFC/',sc_subgroup,'/');
outmat = sFC_cv;
outtab = array2table(outmat);     
outfile = append(sFCpath,'group_cv_sFC.csv');
writetable(outtab,outfile,'WriteRowNames',0,'WriteVariableNames',0) 

%Save group strength for sFC.
outmat = sFC_med;
outtab = array2table(outmat);     
outfile = append(sFCpath,'group_med_sFC.csv');
writetable(outtab,outfile,'WriteRowNames',0,'WriteVariableNames',0) 
end
