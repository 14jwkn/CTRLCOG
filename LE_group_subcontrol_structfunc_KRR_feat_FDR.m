%{
Across participants, for a type of SC normalization, and the type and percentage 
of thresholding, do FDR across the Haufe feature scores from all models of 
interest.
Output:
cctrl,'_',cstatetype,'_',cseptype,'_TwoP_BHFDR_covha_feat.csv' FDR-corrected p-values for the Haufe feature scores.
%}

%Define command line arguments.
function [] = LE_group_subcontrol_structfunc_KRR_feat_FDR(k,sctype,threstype,thresval,controltype,statetype,septype,nrep,inner_k,outer_k,nperm)
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
disp(append('Doing: ',k,' ',sctype,' ',threstype,' ',thresval,' ',controltype,' ',...
            statetype,' ',septype,' ',...
            nrep,' ',inner_k,' ',outer_k,' ',nperm))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set paths.            
subgroup = 'full';
basepath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                  subgroup,'/',k,'/SC_dFC/',...
                  sc_subgroup,'/collect/',threstype,'/',thresval,'/');
fdrpath = append(basepath,'KRRXFS/',controltype,'_',statetype,'_',septype,'_',...
                nrep,'_',inner_k,'_',outer_k,'_',sctype,'/'); 

%Create path.
if not(isfolder(fdrpath))
    mkdir(fdrpath)
end

%Set parameters.
cctrl = 'ave';
cstatetype = 'SC_dFCcat';
cseptype = 'comCFAng';
if strcmp(cseptype,'comCFAng')
    ylabs = string({'gCFA' 'P24_CR' 'PV' 'gFngCFA' 'gCngCFA'});
    ny = size(ylabs,2);
    inpath = append(basepath,'KRRXFS/',cctrl,'_',cstatetype,'_',cseptype,'_',...
                nrep,'_',inner_k,'_',outer_k,'_',sctype,'/');
    infile = append(inpath,'gCFA_covha_featscores.csv');
    cmat = readtable(infile);
    xlabs = string(table2array(cmat(:,1)));
    nx = size(xlabs,1);    
end
ctrlvers = string({'ave','mod'});
nct = size(ctrlvers,2);

%Go through iterations.
pmat = zeros(nct,nx,ny);
for ctidx=1:nct
    cctrl = ctrlvers(ctidx);
    for cidx=1:ny
        ccog = ylabs(cidx);

        %Read in p-values.
        inpath = append(basepath,'KRRXFS/',cctrl,'_',cstatetype,'_',cseptype,'_',...
                        nrep,'_',inner_k,'_',outer_k,'_',sctype,'/');
        infile = append(inpath,ccog,'_covha_featscores.csv');
        inmat = readtable(infile);
        inmat = inmat(:,'TwoP');
        pmat(ctidx,:,cidx) = table2array(inmat);
    end
end

%Convert into one vector.
pvec = reshape(pmat,[(nct*nx*ny) 1]);

%FDR correct for every test.
all_fdrvec = mafdr(pvec,'BHFDR',1);
allfdr = reshape(all_fdrvec,[nct nx ny]);allfdr_thres = allfdr;
allfdr_thres(allfdr >= 0.05) = NaN;

%Save.
for ctidx=1:nct
    cctrl = ctrlvers(ctidx);
    outfile = append(fdrpath,cctrl,'_',cstatetype,'_',cseptype,'_TwoP_BHFDR_covha_feat.csv');
    outtab = array2table(squeeze(allfdr(ctidx,:,:)));
    outtab = [array2table(xlabs) outtab];
    outtab.Properties.VariableNames = ['Label',ylabs];
    writetable(outtab,outfile)
end
disp('Saved.')
end
