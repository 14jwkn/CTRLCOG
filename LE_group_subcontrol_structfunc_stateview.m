%{
Across participants, for a type of SC normalization, and the type and percentage 
of thresholding, average dFC, sFC, and SC matrices across participants and
generate thresholded versions. Save the Gramian eigenvalues, average
controllability, and modal controllability for each matrix. Also save which
regions have errors regarding Gramian eigenvalues, to remove in later
controllability analyses.
Output:
SC_sFC_dFC.h5 Thresholded dFC, sFC, and SC matrices averaged across participants.
dFC_gram.csv Gramian eigenvalues for average dFC across participants.
dFC_ave.csv Average controllability for average dFC across participants.
dFC_mod.csv Modal controllability for average dFC across participants.
sFC_gram.csv Gramian eigenvalues for average sFC across participants.
sFC_ave.csv Average controllability for average sFC across participants.
sFC_mod.csv Modal controllability for average sFC across participants.
sctype,'_SC_gram.csv' Gramian eigenvalues for average SC across participants for a specific normalization.
sctype,'_SC_ave.csv' Average controllability for average SC across participants for a specific normalization.
sctype,'_SC_mod.csv' Modal controllability for average SC across participants for a specific normalization.
sctype,'_SC_sFC_dFC_gram.csv' Records regions with Gramian errors for each average SC, sFC, and dFC matrix.
%}

%Define command line arguments.
function [] = LE_group_subcontrol_structfunc_stateview(k,sctype,threstype,thresval)
% k = Selected k for k-medians clustering.
% sctype = Type of SC normalization.
% threstype = Type of matrix thresholding. 
% thresval = Thresholding percentage.
disp(append('Doing: ',k,sctype,' ',threstype,' ',thresval))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB_packages_only'))

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
pkeep = 1 - (str2num(thresval)/100);
pthrow = str2num(thresval);
expA = ((nroi * nroi) - nroi)/2;

%Read in dFC state values.
infile = append(basepath,'group_med_dFC.csv');
if strcmp(k,'11')
    dFC_state = readmatrix(infile)';
else
    dFC_state = readmatrix(infile);
end

%Read in dFC CV values.         
infile = append(basepath,'group_cv_dFC.csv');
if  strcmp(k,'11')
    dFC_cv = readmatrix(infile)';
else
    dFC_cv = readmatrix(infile);
end

%Read in sFC values.
infile = append('../outputs/r_sFC/',clean,'/',order,'/',roiname,'/',sc_subgroup,'/group_med_sFC.csv');
sFC_state = readmatrix(infile);

%Read in sFC values.
infile = append('../outputs/r_sFC/',clean,'/',order,'/',roiname,'/',sc_subgroup,'/group_cv_sFC.csv');
sFC_cv = readmatrix(infile);

%Read in SC values.
infile = append('../outputs/d_SC/',sc_subgroup,'/group_med_',sctype,'.csv');
sc_state = readmatrix(infile);

%Read in SC CV values.
infile = append('../outputs/d_SC/',sc_subgroup,'/group_cv_',sctype,'.csv');
sc_cv = readmatrix(infile);

%Generate Gramians and controllability for dFC.
grammed = zeros(nk,nroi);
avemed = zeros(nk,nroi);
modmed = zeros(nk,nroi);

%For each centroid for dFC.
for kidx = 1:nk

    %Convert dmed into square.
    disp(append('dFC ',num2str(kidx)))
    sq_dmed = zeros(nroi,nroi);
    sq_dmed(triu(true(nroi),1)) = dFC_state(kidx,:);
    sq_dmed = sq_dmed + sq_dmed.';
    sq_dmed = sq_dmed.*~eye(size(sq_dmed));

    %Threshold by group consistency with minimum spanning tree.
    if strcmp(threstype,'groupconsist_mst')

        %Convert to square.
        sq_mat = zeros(nroi,nroi);
        sq_mat(triu(true(nroi),1)) = dFC_cv(kidx,:);
        sq_cv = sq_mat + sq_mat.';

        %Create MST and iteratively add back the smallest edges until target
        %percentage reached.
        sq_cv_th = thres_with_mst(sq_cv,pkeep);

        %Threshold. 
        th_dmed = sq_dmed;
        th_dmed(sq_cv_th==0) = 0; 
    
    %Threshold by group consistency.
    elseif strcmp(threstype,'groupconsist')
        
        %Generate square CV matrix.
        sq_mat = zeros(nroi,nroi);
        sq_mat(triu(true(nroi),1)) = dFC_cv(kidx,:);
        sq_cv = sq_mat + sq_mat.';
        
        %Threshold.
        th_dmed = threshold_arbmeasure(sq_dmed,-sq_cv,pkeep);
    
    %Threshold by group strength.
    elseif strcmp(threstype,'groupstr')
        
        %Generate square strength matrix.
        sq_mat = zeros(nroi,nroi);
        sq_mat(triu(true(nroi),1)) = dFC_state(kidx,:);
        sq_str = sq_mat + sq_mat.';
        
        %Threshold.
        th_dmed = threshold_arbmeasure(sq_dmed,sq_str,pkeep);
    
    %No threshold.
    elseif strcmp(threstype,'none')
    
        %Set again.
        th_dmed = sq_dmed;
    end    
    
    %Save the matrix.
    outfile = append(outpath,'SC_sFC_dFC.h5');
    cname = append('s',num2str(kidx));
    cvar = th_dmed;
    outkey = append('/',cname); 
    try
        h5create(outfile,outkey,size(cvar));
    catch
        disp('File and key created already.')
    end
    h5write(outfile,outkey,cvar);
    
    %Print network density.
    A = th_dmed;
    triA = nnz(A(triu(true(size(A)),1)));
    disp(append('Density: ',num2str(triA/expA)))
    
    %Go through each node, setting that node as the control node.
    A = th_dmed;
    for ridx = 1:nroi
        
        %Create B matrix.
        B = eye(nroi);
        B = B(:,ridx);

        %Determine minimum Controllability Gramian eigenvalue, 
        %set 0 if not controllable.
        try
            [~,smeig] = Gramian(A,B);
            grammed(kidx,ridx) = smeig;
        catch
            grammed(kidx,ridx) = 0;
        end
    end
    disp(append('Gramian errors: ',num2str(sum(grammed(kidx,:)<=0))))
    
    %Calculate AC and MC.
    A = th_dmed;
    try
        avemed(kidx,:) = averMeas(A);              
    catch
        avemed(kidx,:) = repmat(NaN,360,1);             
    end
    try
        modmed(kidx,:) = moduMeas(A);              
    catch
        modmed(kidx,:) = repmat(NaN,360,1);             
    end    
end
dFC_gram = sum((grammed<=0),2);

%Save Gramian eigenvalues and controllability.
rowlab = append('s',string(1:nk));
ctab = grammed;
outtab = array2table(ctab,'RowNames',rowlab);
outfile = append(outpath,'dFC_gram.csv');
writetable(outtab,outfile,'WriteRowNames',1,'WriteVariableNames',0) 
ctab = avemed;
outtab = array2table(ctab,'RowNames',rowlab);
outfile = append(outpath,'dFC_ave.csv');
writetable(outtab,outfile,'WriteRowNames',1,'WriteVariableNames',0) 
ctab = modmed;
outtab = array2table(ctab,'RowNames',rowlab);
outfile = append(outpath,'dFC_mod.csv');
writetable(outtab,outfile,'WriteRowNames',1,'WriteVariableNames',0) 

%Generate Gramians and controllability for sFC.
disp('sFC')
grammed = zeros(1,nroi);
avemed = zeros(1,nroi);
modmed = zeros(1,nroi);

%Convert dmed into square.
sq_dmed = zeros(nroi,nroi);
sq_dmed(triu(true(nroi),1)) = sFC_state;
sq_dmed = sq_dmed + sq_dmed.';
sq_dmed = sq_dmed.*~eye(size(sq_dmed));

%Threshold by group consistency with minimum spanning tree.
if strcmp(threstype,'groupconsist_mst')

    %Convert to square.
    sq_mat = zeros(nroi,nroi);
    sq_mat(triu(true(nroi),1)) = sFC_cv;
    sq_cv = sq_mat + sq_mat.';

    %Create MST and iteratively add back the smallest edges until target
    %percentage reached.
    sq_cv_th = thres_with_mst(sq_cv,pkeep);

    %Threshold. 
    th_dmed = sq_dmed;
    th_dmed(sq_cv_th==0) = 0;  

%Threshold by group consistency.
elseif strcmp(threstype,'groupconsist')

    %Generate square CV matrix.
    sq_mat = zeros(nroi,nroi);
    sq_mat(triu(true(nroi),1)) = sFC_cv;
    sq_cv = sq_mat + sq_mat.';

    %Threshold.
    th_dmed = threshold_arbmeasure(sq_dmed,-sq_cv,pkeep);

%Threshold by group strength.
elseif strcmp(threstype,'groupstr')

    %Generate square strength matrix.
    sq_mat = zeros(nroi,nroi);
    sq_mat(triu(true(nroi),1)) = sFC_state;
    sq_str = sq_mat + sq_mat.';

    %Threshold.
    th_dmed = threshold_arbmeasure(sq_dmed,sq_str,pkeep);

%No threshold.
elseif strcmp(threstype,'none')

    %Set again.
    th_dmed = sq_dmed;
end    

%Save the matrix.
outfile = append(outpath,'SC_sFC_dFC.h5');
cname = append('sFC');
cvar = th_dmed;
outkey = append('/',cname); 
try
    h5create(outfile,outkey,size(cvar));
catch
    disp('File and key created already.')
end
h5write(outfile,outkey,cvar);

%Print network density.
A = th_dmed;
triA = nnz(A(triu(true(size(A)),1)));
disp(append('Density: ',num2str(triA/expA)))

%Go through each node, setting that node as the control node.
A = th_dmed;
for ridx = 1:nroi

    %Create B matrix.
    B = eye(nroi);
    B = B(:,ridx);

    %Determine minimum Controllability Gramian eigenvalue, 
    %set 0 if not controllable.
    try
        [~,smeig] = Gramian(A,B);
        grammed(1,ridx) = smeig;
    catch
        grammed(1,ridx) = 0;
    end
end
disp(append('Gramian errors: ',num2str(sum(grammed(1,:)<=0))))
sFC_gram = sum((grammed<=0),2);

%Calculate AC and MC.
A = th_dmed;
try
    avemed(1,:) = averMeas(A);              
catch
    avemed(1,:) = repmat(NaN,360,1);             
end
try
    modmed(1,:) = moduMeas(A);              
catch
    modmed(1,:) = repmat(NaN,360,1);             
end 

%Save Gramian eigenvalues and controllability.
rowlab = {'sFC'};
ctab = grammed;
outtab = array2table(ctab,'RowNames',rowlab);
outfile = append(outpath,'sFC_gram.csv');
writetable(outtab,outfile,'WriteRowNames',1,'WriteVariableNames',0) 
ctab = avemed;
outtab = array2table(ctab,'RowNames',rowlab);
outfile = append(outpath,'sFC_ave.csv');
writetable(outtab,outfile,'WriteRowNames',1,'WriteVariableNames',0) 
ctab = modmed;
outtab = array2table(ctab,'RowNames',rowlab);
outfile = append(outpath,'sFC_mod.csv');
writetable(outtab,outfile,'WriteRowNames',1,'WriteVariableNames',0) 

%Generate Gramian for SC.
disp('SC')
grammed = zeros(1,nroi);
avemed = zeros(1,nroi);
modmed = zeros(1,nroi);

%Convert dmed into square.
sq_dmed = zeros(nroi,nroi);
sq_dmed(triu(true(nroi),1)) = sc_state;
sq_dmed = sq_dmed + sq_dmed.';
sq_dmed = sq_dmed.*~eye(size(sq_dmed));

%Threshold by group consistency with minimum spanning tree.
if strcmp(threstype,'groupconsist_mst')

    %Convert to square.
    sq_mat = zeros(nroi,nroi);
    sq_mat(triu(true(nroi),1)) = sc_cv;
    sq_cv = sq_mat + sq_mat.';

    %Create MST and iteratively add back the smallest edges until target
    %percentage reached.
    sq_cv_th = thres_with_mst(sq_cv,pkeep);

    %Threshold. 
    th_dmed = sq_dmed;
    th_dmed(sq_cv_th==0) = 0;  

%Threshold by group consistency.
elseif strcmp(threstype,'groupconsist')

    %Generate square CV matrix.
    sq_mat = zeros(nroi,nroi);
    sq_mat(triu(true(nroi),1)) = sc_cv;
    sq_cv = sq_mat + sq_mat.';

    %Threshold.
    th_dmed = threshold_arbmeasure(sq_dmed,-sq_cv,pkeep);

%Threshold by group strength.
elseif strcmp(threstype,'groupstr')

    %Generate square strength matrix.
    sq_mat = zeros(nroi,nroi);
    sq_mat(triu(true(nroi),1)) = sc_state;
    sq_str = sq_mat + sq_mat.';

    %Threshold.
    th_dmed = threshold_arbmeasure(sq_dmed,sq_str,pkeep);

%No threshold.
elseif strcmp(threstype,'none')

    %Set again.
    th_dmed = sq_dmed;
end    

%Save the matrix.
outfile = append(outpath,'SC_sFC_dFC.h5');
cname = append('SC_',sctype);
cvar = th_dmed;
outkey = append('/',cname); 
try
    h5create(outfile,outkey,size(cvar));
catch
    disp('File and key created already.')
end
h5write(outfile,outkey,cvar);

%Print network density.
A = th_dmed;
triA = nnz(A(triu(true(size(A)),1)));
disp(append('Density: ',num2str(triA/expA)))

%Go through each node, setting that node as the control node.
A = th_dmed;
for ridx = 1:nroi

    %Create B matrix.
    B = eye(nroi);
    B = B(:,ridx);

    %Determine minimum Controllability Gramian eigenvalue, 
    %set 0 if not controllable.
    try
        [~,smeig] = Gramian(A,B);
        grammed(1,ridx) = smeig;
    catch
        grammed(1,ridx) = 0;
    end
    
    %Calculate AC and MC.
    A = th_dmed;
    try
        avemed(1,:) = averMeas(A);              
    catch
        avemed(1,:) = repmat(NaN,360,1);             
    end
    try
        modmed(1,:) = moduMeas(A);              
    catch
        modmed(1,:) = repmat(NaN,360,1);             
    end 
end
disp(append('Gramian errors: ',num2str(sum(grammed(1,:)<=0))))
sc_gram = sum((grammed<=0),2);

%Save Gramian eigenvalues and controllability.
rowlab = {'SC'};
ctab = grammed;
outtab = array2table(ctab,'RowNames',rowlab);
outfile = append(outpath,sctype,'_SC_gram.csv');
writetable(outtab,outfile,'WriteRowNames',1,'WriteVariableNames',0) 
ctab = avemed;
outtab = array2table(ctab,'RowNames',rowlab);
outfile = append(outpath,sctype,'_SC_ave.csv');
writetable(outtab,outfile,'WriteRowNames',1,'WriteVariableNames',0) 
ctab = modmed;
outtab = array2table(ctab,'RowNames',rowlab);
outfile = append(outpath,sctype,'_SC_mod.csv');
writetable(outtab,outfile,'WriteRowNames',1,'WriteVariableNames',0) 

%Collect Gramian errors for SC, sFC, and dFC and save.
rowlab = string([{'SC','sFC'} cellstr(append('s',string(1:nk)))]);
ctab = [sc_gram;sFC_gram;dFC_gram];
outtab = array2table(ctab,'RowNames',rowlab);
outfile = append(outpath,sctype,'_SC_sFC_dFC_gram.csv');
writetable(outtab,outfile,'WriteRowNames',1,'WriteVariableNames',0) 
end
