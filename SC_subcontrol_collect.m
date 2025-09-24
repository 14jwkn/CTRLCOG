%{
Across participants, for a type of SC normalization, and the type and percentage of 
thresholding, collect the degree, average controllability, and modal controllability
for SC.
Output:
cver,'_sc.csv' Degree, average controllability, or modal controllability across participants.
%}

%Define command line arguments.
function [] = SC_subcontrol_collect(sctype,threstype,thresval)
% sctype = Type of SC normalization.
% threstype = Type of matrix thresholding. 
% thresval = Thresholding percentage.
disp(append('Doing: ',sctype,' ',threstype,' ',thresval))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set path.
sc_subgroup = 'dr_full';
dtipath = append('../outputs/d_SC/',sc_subgroup,'/');

%Read in subjects.
if strcmp(sc_subgroup,'dr_full') 
    subfile = 'dr_full_intersect.txt';
end        
subjects = string(textread(subfile,'%s','delimiter',','));
nsub = size(subjects,1);

%Initialize values.
nroi = 360;

%Go through each version.
verlist = string({'deg','orig_ave','orig_mod'});
nver = size(verlist,2);
for vidx=1:nver
    
    %Extract.
    cver = verlist(vidx);
    inkey = append('/',cver);
    
    %Go through subjects.
    vermat = zeros(nsub,nroi);
    for sidx=1:nsub
        
        %Extract.
        csub = subjects(sidx);
        inpath = append(dtipath,threstype,'/',csub,'/',thresval,'/',sctype,'/');
        infile = append(inpath,'sub_sc.h5');
        inmat = h5read(infile,inkey);
        
        %Append.
        vermat(sidx,:) = inmat';
    end    
    
    %Save.
    outpath = append(dtipath,threstype,'/',thresval,'/',sctype,'/');
    if not(isfolder(outpath))
        mkdir(outpath)
    end
    outfile = append(outpath,cver,'_sc.csv');
    outtab = array2table(vermat,'RowNames',subjects);     
    writetable(outtab,outfile,'WriteRowNames',1,'WriteVariableNames',0)
end
disp('Saved.')
end
