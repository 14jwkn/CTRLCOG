%{
For each participant, for a k for clustering, and the type and percentage of 
thresholding, generate the thresholded dFC state matrices, degree, average
controllability, and modal controllability.
Output:
'control_thres_med_dFC.h5' h5 file with thresholded subject degree, average controllability, and modal controllability.
%}

%Define command line arguments.
function [] = LE_group_subcontrol_structfunc(k,subject,threstype,thresval)
% k = Selected k for k-medians clustering.
% subject = Subject ID.
% threstype = Type of matrix thresholding. 
% thresval = Thresholding percentage.
disp(append('Doing: ',k,' ',subject,' ',threstype,' ',thresval))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O.       
subgroup = 'full';
sc_subgroup = 'dr_full';
statepath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                   subgroup,'/',k,'/SC_dFC/',sc_subgroup,'/'); 
outpath = append(statepath,subject,'/',threstype,'/',thresval,'/');           
thresmed_file = append(outpath,'control_thres_med_dFC.h5');

%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end

%Initialize values.
nk = str2num(k);
nroi = 360;
pkeep = 1 - (str2num(thresval)/100);
pthrow = str2num(thresval);
expA = ((nroi * nroi) - nroi)/2;

%Get centroid resorting.
classfile = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                       'visclass/',subgroup,'/allk/',subgroup,'_sortclass.csv');
classmat = readtable(classfile,'ReadRowNames',true,'ReadVariableNames',true);
classmat.Properties.VariableNames = string(1:size(classmat,2));
classmat = classmat(k,:);
classmat = rmmissing(classmat,2);

%Read in subject centroids and resort.
inpath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                subgroup,'/',k,'/',subject,'/');
infile = append(inpath,'med_dFC.h5'); 
inkey = append('/',subgroup);
dFCmed = h5read(infile,inkey);
dFCmed = dFCmed(table2array(classmat),:);
        
%Read in group CV values.         
infile = append(statepath,'group_cv_dFC.csv');
if strcmp(k,'11')
    kcv_med = readmatrix(infile)';
else
    kcv_med = readmatrix(infile);
end

%Read in group strength values.
infile = append(statepath,'group_med_dFC.csv');
if strcmp(k,'11')
    kstr_med = readmatrix(infile)';
else
    kstr_med = readmatrix(infile);
end

%Generate AC.
avemed = zeros(nk,nroi);

%Generate MC.
modmed = zeros(nk,nroi);

%Generate degree.
degmed = zeros(nk,nroi);

%For each centroid.
for kidx = 1:nk

    %Convert dmed into square.
    sq_dmed = zeros(nroi,nroi);
    sq_dmed(triu(true(nroi),1)) = dFCmed(kidx,:);
    sq_dmed = sq_dmed + sq_dmed.';
    sq_dmed = sq_dmed.*~eye(size(sq_dmed));
    
    %Threshold by group consistency.
    if strcmp(threstype,'groupconsist')
        
        %Generate square CV matrix.
        sq_mat = zeros(nroi,nroi);
        sq_mat(triu(true(nroi),1)) = kcv_med(kidx,:);
        sq_cv = sq_mat + sq_mat.';
        
        %Threshold.
        th_dmed = threshold_arbmeasure(sq_dmed,-sq_cv,pkeep);
    
    %Threshold by group strength.
    elseif strcmp(threstype,'groupstr')
        
        %Generate square strength matrix.
        sq_mat = zeros(nroi,nroi);
        sq_mat(triu(true(nroi),1)) = kstr_med(kidx,:);
        sq_str = sq_mat + sq_mat.';
        
        %Threshold.
        th_dmed = threshold_arbmeasure(sq_dmed,sq_str,pkeep);
    
    %No threshold.
    elseif strcmp(threstype,'none')
    
        %Set again.
        th_dmed = sq_dmed;
    end    
    
    %Print network density.
    A = th_dmed;
    triA = nnz(A(triu(true(size(A)),1)));
    disp(append('Density: ',num2str(triA/expA)))

    %Get absolute degree.
    degmed(kidx,:) = sum(abs(th_dmed)); 
    
    %Get average controllability.
    try
        avemed(kidx,:) = averMeas(th_dmed);              
    catch
        avemed(kidx,:) = repmat(NaN,360,1);             
    end
    
    %Get modal controllability.
    try
        modmed(kidx,:) = moduMeas(th_dmed);
    catch
        modmed(kidx,:) = repmat(NaN,360,1);
    end
end    

%Generate list to save. 
varlist = {avemed,modmed,degmed};
varnames = string({'ave','mod','deg'});                            
nvars = size(varnames,2);
outfile = thresmed_file;
for vidx=1:nvars
    cname = varnames{vidx};
    cvar = varlist{vidx};
   
    outkey = append('/',cname); 
    try
        h5create(outfile,outkey,size(cvar));
    catch
        disp('File and key created already.')
    end
    h5write(outfile,outkey,cvar);
end
disp('Saved.')
end
