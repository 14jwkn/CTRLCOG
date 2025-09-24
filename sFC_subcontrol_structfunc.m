%{
For each participant, for the type and percentage of thresholding, generate 
the thresholded sFC matrix, degree, average controllability, and modal controllability.
Output:
sub_sFC.h5 h5 file with thresholded subject sFC matrix, degree, average controllability, and modal controllability.
%}

%Define command line arguments.
function [] = sFC_subcontrol_structfunc(subject,threstype,thresval)
% subject = Subject ID.
% threstype = Type of matrix thresholding. 
% thresval = Thresholding percentage.
disp(append('Doing: ',subject,' ',threstype,' ',thresval))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O.       
sc_subgroup = 'dr_full';
sFCpath = append('../outputs/r_sFC/',sc_subgroup,'/');
outpath = append(sFCpath,subject,'/',threstype,'/',thresval,'/');

%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end

%Read in subjects for index.
subfile = 'dr_full_intersect.txt';      
subjects = textread(subfile,'%s','delimiter',',');
subjects = string(subjects); 
subidx = find(strcmp(subjects,subject));

%Initialize values.
nroi = 360;
pkeep = 1 - (str2num(thresval)/100);
pthrow = str2num(thresval);
expA = ((nroi * nroi) - nroi)/2;

%Read in subject sFC.
infile = append(sFCpath,'sub_sFC.h5');
inkey = append('/sFC');
sFCmat = h5read(infile,inkey);
sFCmat = sFCmat(subidx,:);
        
%Read in group CV values.         
infile = append(sFCpath,'group_cv_sFC.csv');
cv_med = readmatrix(infile);

%Read in group strength values.
infile = append(sFCpath,'group_med_sFC.csv');
str_med = readmatrix(infile);

%Convert dmed into square.
sq_dmed = zeros(nroi,nroi);
sq_dmed(triu(true(nroi),1)) = sFCmat;
sq_dmed = sq_dmed + sq_dmed.';
sq_dmed = sq_dmed.*~eye(size(sq_dmed));

%Threshold by group consistency.
if strcmp(threstype,'groupconsist')

    %Generate square CV matrix.
    sq_mat = zeros(nroi,nroi);
    sq_mat(triu(true(nroi),1)) = cv_med;
    sq_cv = sq_mat + sq_mat.';

    %Threshold.
    th_dmed = threshold_arbmeasure(sq_dmed,-sq_cv,pkeep);

%Threshold by group strength.
elseif strcmp(threstype,'groupstr')

    %Generate square strength matrix.
    sq_mat = zeros(nroi,nroi);
    sq_mat(triu(true(nroi),1)) = str_med;
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
degmed = sum(abs(th_dmed))'; 

%Get average controllability.
try
    avemed = averMeas(th_dmed);              
catch
    avemed = repmat(NaN,360,1);             
end

%Get modal controllability.
try
    modmed = moduMeas(th_dmed);
catch
    modmed = repmat(NaN,360,1);
end

%Generate list to save.
varlist = {avemed,modmed,degmed};
varnames = string({'ave','mod','deg'});         
nvars = size(varnames,2);
outfile = append(outpath,'sub_sFC.h5');
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
