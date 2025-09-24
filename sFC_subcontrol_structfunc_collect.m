%{
Across participants, for the type and percentage of thresholding, collect the 
degree, average controllability, and modal controllability for sFC.
Output:
cver,'_sFC.csv' Degree, average controllability, or modal controllability across participants.
%}

%Define command line arguments.
function [] = sFC_subcontrol_structfunc_collect(threstype,thresval)
% threstype = Type of matrix thresholding. 
% thresval = Thresholding percentage.
disp(append('Doing: ',threstype,' ',thresval))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O.       
sc_subgroup = 'dr_full';
sFCpath = append('../outputs/r_sFC/',sc_subgroup,'/');
outpath = append(sFCpath,'/',threstype,'/',thresval,'/');  
if not(isfolder(outpath))
    mkdir(outpath)
end

%Read in subjects.
subfile = 'dr_full_intersect.txt';      
subjects = string(textread(subfile,'%s','delimiter',','));
nsub = size(subjects,1);

%Initialize values.
nroi = 360;

%Go through each version.
verlist = string({'ave','mod','deg'});
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
        inpath = append(sFCpath,'/',csub,'/',threstype,'/',thresval,'/');
        infile = append(inpath,'sub_sFC.h5');
        inmat = h5read(infile,inkey);
        
        %Append.
        vermat(sidx,:) = inmat';
    end    
    
    %Save.
    outfile = append(outpath,cver,'_sFC.csv');
    outtab = array2table(vermat,'RowNames',subjects);     
    writetable(outtab,outfile,'WriteRowNames',1,'WriteVariableNames',0)
end
end
