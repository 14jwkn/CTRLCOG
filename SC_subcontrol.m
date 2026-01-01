%{
For each participant, for a type of SC normalization, and the type and percentage of 
thresholding, generate the thresholded SC matrix, degree, average
controllability, and modal controllability.
Output:
'sub_sc.h5' h5 file with thresholded subject SC matrix, degree, average controllability, and modal controllability.
%}

%Define command line arguments.
function [] = SC_subcontrol(sctype,subject,threstype,thresval)
% sctype = Type of SC normalization.
% subject = Subject ID.
% threstype = Type of matrix thresholding. 
% thresval = Thresholding percentage.
disp(append('Doing: ',sctype,' ',subject,' ',threstype,' ',thresval))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set path.
sc_subgroup = 'dr_full';
basepath = append('../inputs/data/dti/');
dtipath = append('../outputs/d_SC/',sc_subgroup,'/');

%Initialize values.
thresnum = str2num(thresval);
pkeep = 1 - (thresnum/100);
pthrow = thresnum;
nroi = 360;

%Read in subject matrix.
inkey = append('/',subject);  
subpath = append(basepath,subject,'/');
infile = append(subpath,'sub_sc_',sctype,'.h5');
submat = h5read(infile,inkey); 

%Remove diagonal.
submat = submat.*~eye(size(submat));

%If doing group CV thresholding with minimum spanning tree.
if strcmp(threstype,'groupconsist_mst')

    %Read in group CV matrix.
    infile = append(dtipath,'group_cv_',sctype,'.csv');
    groupcv = readmatrix(infile);
    
    %Convert to square.
    sq_cv = zeros(nroi,nroi);
    sq_cv(triu(true(nroi),1)) = groupcv;
    sq_cv = sq_cv + sq_cv.';
    
    %Create MST and iteratively add back the smallest edges until target
    %percentage reached.
    sq_cv_th = thres_with_mst(sq_cv,pkeep);
    
    %Threshold subject matrix with it.      
    submat(sq_cv_th==0) = 0;  

%If doing group CV thresholding.
elseif strcmp(threstype,'groupconsist')
    
    %Read in group CV matrix.
    infile = append(dtipath,'group_cv_',sctype,'.csv');
    groupcv = readmatrix(infile);
    
    %Convert to square.
    sq_cv = zeros(nroi,nroi);
    sq_cv(triu(true(nroi),1)) = groupcv;
    sq_cv = sq_cv + sq_cv.';
    
    %Threshold subject matrix with it.
    submat = threshold_arbmeasure(submat,-sq_cv,pkeep);

%If doing group strength thresholding.
elseif strcmp(threstype,'groupstr')
    
    %Read in group strength matrix.
    infile = append(dtipath,'group_med_',sctype,'.csv');
    groupcv = readmatrix(infile);
    
    %Convert to square.
    sq_cv = zeros(nroi,nroi);
    sq_cv(triu(true(nroi),1)) = groupcv;
    sq_cv = sq_cv + sq_cv.';
    
    %Threshold subject matrix with it.
    submat = threshold_arbmeasure(submat,sq_cv,pkeep);

%If doing subject strength thresholding.
elseif strcmp(threstype,'substr')
    submat = threshold_strength(submat,pthrow);
end

%Save matrix.
outkey = append('/',subject); 
outpath = append(dtipath,threstype,'/',subject,'/',thresval,'/',sctype,'/');
if not(isfolder(outpath))
    mkdir(outpath)
end
outfile = append(outpath,'sub_sc.h5');
try
    h5create(outfile,outkey,size(submat));
catch
    disp('File and key created already.')
end
h5write(outfile,outkey,submat);

%Calculate absolute degree.
deg = sum(abs(submat)); 

%Calculate AC/MC with discrete model normalization from original paper.
try
    orig_ave = averMeas(submat);
catch
    orig_ave = repmat(NaN,360,1);
end
try
    orig_mod = moduMeas(submat);
catch
    orig_mod = repmat(NaN,360,1);
end

%Save these values.
outlist = string({'deg','orig_ave','orig_mod'});
nout = size(outlist,2);
outfile = append(outpath,'sub_sc.h5');
for oidx=1:nout
    cout = outlist(oidx);
    outkey = append('/',cout); 
    if strcmp(cout,'deg')
        cvar = deg;
    elseif strcmp(cout,'orig_ave') 
        cvar = orig_ave;
    elseif strcmp(cout,'orig_mod') 
        cvar = orig_mod;
    end
    try
        h5create(outfile,outkey,size(cvar));
    catch
        disp('File and key created already.')
    end
    h5write(outfile,outkey,cvar);
end
end
