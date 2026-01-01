%{
Across participants, for a type of SC normalization, and the type and percentage 
of thresholding, do FDR across the test set CV accuracy scores and accuracy
comparison scores of interest.
Output:
testacc_corrected.csv FDR-corrected p-values for test set CV accuracy scores.
ccompset,'_corrected.csv' FDR-corrected p-values for comparison scores.
%}

%Define command line arguments.
function [] = LE_group_subcontrol_structfunc_KRR_scorcomp_FDR(k,sctype,threstype,thresval,controltype,statetype,septype,nrep,inner_k,outer_k,nperm)
% k = Selected k for k-medians clustering.
% sctype = Type of SC normalization.
% threstype = Type of matrix thresholding. 
% thresval = Thresholding percentage.
% controltype = Control type.
% statetype = Which states.
% septype = Cognitive variables batch.
% nrep = Number of repetitions for the CV.
% inner_k = Number of inner loops for the CV.
% outer_k = Number of outer loops for the CV.
% nperm = Number of permutations for permutation testing.
disp(append('Doing: ',k,' ',...
            sctype,' ',threstype,' ',thresval,' ',controltype,' ',...
            statetype,' ',septype,' ',...
            nrep,' ',inner_k,' ',outer_k))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set paths and make missing ones.         
subgroup = 'full';
sc_subgroup = 'dr_full';
basepath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                  subgroup,'/',k,'/SC_dFC/',...
                  sc_subgroup,'/collect/',threstype,'/',thresval,'/');
fpath = append(basepath,'KRRXF/',controltype,'_',statetype,'_',septype,'_',...
               nrep,'_',...
               inner_k,'_',outer_k,'_',sctype,'/');  
fspath = append(basepath,'KRRXFS/',controltype,'_',statetype,'_',septype,'_',...
                nrep,'_',...
                inner_k,'_',outer_k,'_',sctype,'/');      
if not(isfolder(fspath))
    mkdir(fspath)
end           

%Set parameters to read in the overall prediction p-values.
cctrl_list =  string({'ave' 'mod'});
nctrl = size(cctrl_list,2);
cstate_list = string({'SC' 'dFCcat' 'SC_dFCcat'});
nstate = size(cstate_list,2);
cseptype = 'comCFAng';
ccoglist = string({'gCFA' 'P24_CR' 'PV' 'gFngCFA' 'gCngCFA'});
ncog = size(ccoglist,2);

%Go through each matrix and place it in.
predp = zeros(nctrl*ncog*nstate,1);
startidx = 1;
for cctrl=cctrl_list
    for cstate=cstate_list
        
        %Read.
        inpath = append(basepath,'KRRXFS/',cctrl,'_',cstate,'_',cseptype,'_',...
                        nrep,'_',...
                        inner_k,'_',outer_k,'_',sctype,'/');
        infile = append(inpath,'testacc.csv');
        inmat = readtable(infile,ReadRowNames=true);
        
        %Increment.
        endidx = startidx + size(inmat,1) - 1;
        
        %Place it in.
        predp(startidx:endidx,1) = inmat{:,'OneP'};
        
        %Increment.
        startidx = endidx + 1;
    end
end

%Set parameters to read in the overall prediction comparison p-values.
compsets = string({'dFC_vs_SC' 'SC_dFC_vs_dFC' 'cogcompare' 'strcompare'});
ncompsets = size(compsets,2);

%Retrieve.
compp = zeros(0,1);
for coidx=1:ncompsets
    
    %Read.
    ccompset = compsets(coidx);
    infile = append(fpath,ccompset,'.csv');
    inmat = readtable(infile,ReadRowNames=true);
    cpvec = inmat{:,'P'};
    
    %Place it in.
    compp = vertcat(compp,cpvec);
end

%Generate full p vector.
pvec = vertcat(predp,compp);

%FDR correct for every test.
bhvec = mafdr(pvec,'BHFDR',1);

%Go through each matrix and slice out.
startidx = 1;
for cctrl=cctrl_list
    for cstate=cstate_list
        
        %Read.
        inpath = append(basepath,'KRRXFS/',cctrl,'_',cstate,'_',cseptype,'_',...
                        nrep,'_',inner_k,'_',outer_k,'_',sctype,'/');
        infile = append(inpath,'testacc.csv');
        inmat = readtable(infile,ReadRowNames=true);
        
        %Increment.
        endidx = startidx + size(inmat,1) - 1;
        
        %Extract the FDR corrected p and place it in.
        cfdr_vec = bhvec;
        cfdr_lab = 'BH';
        cpart = cfdr_vec(startidx:endidx);
        inmat(:,cfdr_lab) = array2table(cpart);
        
        %Save it.
        outfile = append(inpath,'testacc_corrected.csv');
        writetable(inmat,outfile,WriteRowNames=true);
        
        %Increment.
        startidx = endidx + 1;
    end
end

%Go through each matrix and slice out.
for coidx=1:ncompsets
    
    %Read.
    ccompset = compsets(coidx);
    infile = append(fpath,ccompset,'.csv');
    inmat = readtable(infile,ReadRowNames=true);
    
    %Increment.
    endidx = startidx + size(inmat,1) - 1;
    
    %Extract the FDR corrected p and place it in.
    cfdr_vec = bhvec;
    cfdr_lab = 'BH';
    cpart = cfdr_vec(startidx:endidx);
    inmat(:,cfdr_lab) = array2table(cpart);

    %Save it.
    outfile = append(fspath,ccompset,'_corrected.csv');
    writetable(inmat,outfile,WriteRowNames=true);

    %Increment.
    startidx = endidx + 1;
end
disp('Saved.')
end
