%{
For a k for clustering, and the type and percentage of thresholding, for
each thresholded dFC state collect degree, average controllability, and modal 
controllability across participants.
Output:
cname,'_tab.csv' Degree, average controllability, and modal controllability across participants.
%}

%Define command line arguments.
function [] = LE_group_subcontrol_structfunc_collect(k,threstype,thresval)
% k = Selected k for k-medians clustering.
% threstype = Type of matrix thresholding. 
% thresval = Thresholding percentage.
disp(append('Doing: ',k,' ',threstype,' ',thresval))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up paths.          
subgroup = 'full';
sc_subgroup = 'dr_full';
statepath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                   subgroup,'/',k,'/SC_dFC/',sc_subgroup,'/'); 
outpath = append(statepath,'collect/',threstype,'/',thresval,'/');   

%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end

%Read in subjects.
subfile = 'dr_full_intersect.txt';
subjects = textread(subfile,'%s','delimiter',',');
subjects = string(subjects); 
nsubj = size(subjects,1);

%Set parameters.
nk = str2num(k);
nroi = 360;  
varnames = string({'ave','mod','deg'});                 
nvars = size(varnames,2);

%Go through each matrix type.
for vidx=1:nvars
    
    %Set table.
    ctab = zeros(nsubj,nroi*nk);
    
    %Extract.
    cname = varnames{vidx};
    disp(cname)
    
    %Go through states.
    for kidx = 1:nk
        
        %Set up indices.
        startk = (kidx*nroi) - nroi + 1;
        endk = startk + nroi - 1;
        
        %For each subject.
        for sidx = 1:nsubj

            %Extract.
            csub = subjects(sidx);

            %Set in files.
            inpath = append(statepath,csub,'/',threstype,'/',thresval,'/');   
            infile = append(inpath,'control_thres_med_dFC.h5');

            %Read in and append.
            inkey = append('/',cname);
            inmat = h5read(infile,inkey);
            ctab(sidx,startk:endk) = inmat(kidx,:);
        end
    end  
    
    %Save the table.
    outfile = append(outpath,cname,'_tab.csv');
    outtab = array2table(ctab,'RowNames',subjects);     
    writetable(outtab,outfile,'WriteRowNames',1,'WriteVariableNames',0)  
end
disp('Saved.')
end
