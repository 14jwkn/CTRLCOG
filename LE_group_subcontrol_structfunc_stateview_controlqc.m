%{
For a type of SC normalization, and the type and percentage 
of thresholding, for the average dFC, sFC, and SC matrices across
participants, find the regions with errors regarding Gramian eigenvalues, 
verify the error in their corresponding values by collecting the Gramian 
eigenvalue, average controllability, and modal controllability. 
Output:
cstate,'_err.csv' Gramian eigenvalue, average controllability, and modal controllability for regions with errors regarding Gramian eigenvalues.
%}

%Define command line arguments.
function [] = LE_group_subcontrol_structfunc_stateview_controlqc(k,sctype,threstype,thresval)
% k = Selected k for k-medians clustering.
% sctype = Type of SC normalization.
% threstype = Type of matrix thresholding. 
% thresval = Thresholding percentage.
disp(append('Doing: ',k,sctype,' ',threstype,' ',thresval))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set up I/O.       
subgroup = 'full';
sc_subgroup = 'dr_full';
basepath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                  subgroup,'/',k,'/SC_dFC/',sc_subgroup,'/'); 
outpath = append(basepath,'collect/',threstype,'/',thresval,'/state_images/'); 
    
%If the output folder is not created, create it.
if not(isfolder(outpath))
    mkdir(outpath)
end

%Initialize values.
nk = str2num(k);
nroi = 360;
statelist = string([{append(sctype,'_SC'),'sFC'} cellstr(append('s',string(1:nk)))]);
nstate = size(statelist,2);

%Read in values.
grammat = zeros(2+nk,361);
infile = append(outpath,sctype,'_SC_gram.csv');
grammat(1,:) = readmatrix(infile);
infile = append(outpath,'sFC_gram.csv');
grammat(2,:) = readmatrix(infile);
infile = append(outpath,'dFC_gram.csv');
grammat(3:(nk+2),:) = readmatrix(infile);
grammat = grammat(:,2:end);
avemat = zeros(2+nk,361);
infile = append(outpath,sctype,'_SC_ave.csv');
avemat(1,:) = readmatrix(infile);
infile = append(outpath,'sFC_ave.csv');
avemat(2,:) = readmatrix(infile);
infile = append(outpath,'dFC_ave.csv');
avemat(3:(nk+2),:) = readmatrix(infile);
avemat = avemat(:,2:end);
modmat = zeros(2+nk,361);
infile = append(outpath,sctype,'_SC_mod.csv');
modmat(1,:) = readmatrix(infile);
infile = append(outpath,'sFC_mod.csv');
modmat(2,:) = readmatrix(infile);
infile = append(outpath,'dFC_mod.csv');
modmat(3:(nk+2),:) = readmatrix(infile);
modmat = modmat(:,2:end);

%For each state, find the regions with errors and their corresponding
%control values.
for stidx=1:nstate
    
    %Extract.
    cstate = statelist(stidx);
    cgram = grammat(stidx,:);
    cave = avemat(stidx,:);
    cmod = modmat(stidx,:);
    
    %Find regions.
    reg_err = find(cgram<=0);
    gram_err = cgram(reg_err);
    ave_err = cave(reg_err);
    mod_err = cmod(reg_err);
    
    %Stack and save if errors exist.
    if ~isempty(reg_err)
        tot_err = [reg_err;gram_err;ave_err;mod_err]';
        collab = string({'Region','Gramian','AC','MC'});
        outtab = array2table(tot_err,'VariableNames',collab);
        outfile = append(outpath,cstate,'_err.csv');
        writetable(outtab,outfile,'WriteRowNames',0,'WriteVariableNames',1) 
    end
end 
end
