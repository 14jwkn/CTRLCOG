%{
Across participants, for a type of SC normalization, and the type and percentage 
of thresholding, retrieve the brain map analysis statistically tested values 
of interest to FDR correct. These include the Dice coefficients between 
significant Haufe scores and PFIT/MD maps and R2 for quadratic regression to 
quantify the statistical relationship between Haufe scores and 
controllability/principal cortical gradient. 
Output:
clab,'.csv' Gathered values, p-values, FDR-corrected p-values, and thresholded values based on significance.
%}

%Define command line arguments.
function [] = LE_group_subcontrol_structfunc_KRR_maplook_FDR(k,sctype,threstype,thresval,controltype,statetype,septype,nrep,inner_k,outer_k,wantperm)
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
% wantperm = What number of permutations to run.
disp(append('Doing: ',k,' ',sctype,' ',threstype,' ',thresval,' ',controltype,' ',...
            statetype,' ',septype,' ',...
            nrep,' ',inner_k,' ',outer_k,' ',...
            wantperm))

%Add personal folder to the MATLAB path.
addpath(genpath('../MATLAB'))

%Set paths.
subgroup = 'full';
sc_subgroup = 'dr_full';
basepath = append('../outputs/r_stateflex/statecalc_test/LE/ver_MATLAB/group/',...
                  subgroup,'/',k,'/SC_dFC/',...
                  sc_subgroup,'/collect/',threstype,'/',thresval,'/');
explorepath = append(basepath,'KRRXFS/',controltype,'_',statetype,'_',...
                     nrep,'_',inner_k,'_',outer_k,'_',sctype,'/correxplore/specialp/');     
                 
%Set parameters to read.
nk = str2num(k);
ctrl_list =  string({'ave' 'mod'});
nctrl = size(ctrl_list,2);
block_list = string({'areaDC','ctrlpred','grad1pred'});
if strcmp(cstatetype,'SC_dFCcat')
    statelabs = ['sc' append('s',string(1:nk))];
    ncols = size(statelabs,2);
end
corrtype = 'quadreg';
outpath = append(explorepath,'/FDR_',corrtype,'_',wantperm,'/');
if not(isfolder(outpath))
   mkdir(outpath)
end                                       

%Generate a full matrix.
nrows = (2*2) + 4 + 2;
fulls = zeros(nrows,ncols);
fullp = zeros(nrows,ncols);
full_lab = {};

%Place in the values.
cidx = 1;
for cblock=block_list
    
    %Set input path.
    inpath = append(explorepath,cblock,'/');
    
    %Do areaDC.
    if strcmp(cblock,'areaDC')
        for cctrl=ctrl_list
            
             %Read.
             ccog = 'gCFA';
             infile = append(inpath,'dc_',cctrl,'_',statetype,'_',ccog,'_covha_round.csv');
             inmat = readtable(infile);
             inmat.Properties.RowNames = string(inmat{:,1}); 
             inmat(:,1) = [];
             inmat = inmat(:,{'PFIT','MD'});
             infilep = append(inpath,'dc_',cctrl,'_',statetype,'_',ccog,'_covha_pval.csv');
             inmatp = readtable(infilep);
             inmatp.Properties.RowNames = string(inmatp{:,1}); 
             inmatp(:,1) = [];
             inmatp = inmatp(:,{'PFIT','MD'});
             
             %Place in. 
             fulls(cidx,:) = table2array(inmat(statelabs,'PFIT'));
             fullp(cidx,:) = table2array(inmatp(statelabs,'PFIT'));
             full_lab{cidx} = append('PFIT_',cctrl,'_',ccog);
             cidx = cidx + 1;
             fulls(cidx,:) = table2array(inmat(statelabs,'MD'));
             fullp(cidx,:) = table2array(inmatp(statelabs,'MD'));
             full_lab{cidx} = append('MD_',cctrl,'_',ccog);
             cidx = cidx + 1;
        end  
    
    %Do ctrlpred.
    elseif strcmp(cblock,'ctrlpred')
        cog_list = string({'PV','P24_CR','gCFA'});
        for ccog=cog_list
            if strcmp(ccog,'PV')
                c_ctrl_list = string({'ave'});
            elseif strcmp(ccog,'P24_CR')
                c_ctrl_list = string({'mod'});
            elseif strcmp(ccog,'gCFA')
                c_ctrl_list = string({'ave','mod'});
            end
            for cctrl=c_ctrl_list
                 infile = append(inpath,wantpred,'/',wantpred,'_',cctrl,'_',cstatetype,'_',cseptype,'_',corrtype,'_covha_',cblock,'_round.csv');
                 inmat = readtable(infile);
                 inmat.Properties.RowNames = string(inmat{:,1}); 
                 inmat(:,1) = [];
                 inmat = inmat(ccog,:);
                 infilep = append(inpath,wantpred,'/',wantperm,'_',wantpred,'_',cctrl,'_',cstatetype,'_',cseptype,'_',corrtype,'_covha_',cblock,'_pval.csv');
                 inmatp = readtable(infilep);
                 inmatp.Properties.RowNames = string(inmatp{:,1}); 
                 inmatp(:,1) = [];
                 inmatp = inmatp(ccog,:);
                 fulls(cidx,:) = table2array(inmat);
                 fullp(cidx,:) = table2array(inmatp);
                 full_lab{cidx} = append(cblock,'_',cctrl,'_',ccog);
                 cidx = cidx + 1;
            end 
        end
    
    %Do grad1pred.
    elseif strcmp(cblock,'grad1pred')
        ccog = 'gCFA';
        for cctrl=ctrl_list
             infile = append(inpath,wantpred,'/',cctrl,'_',statetype,'_',septype,'_',corrtype,'_covha_',cblock,'_round.csv');
             inmat = readtable(infile);
             inmat.Properties.RowNames = string(inmat{:,1}); 
             inmat(:,1) = [];
             inmat = inmat(ccog,:);
             infilep = append(inpath,wantpred,'/',wantperm,'_',cctrl,'_',statetype,'_',septype,'_',corrtype,'_covha_',cblock,'_pval.csv');
             inmatp = readtable(infilep);
             inmatp.Properties.RowNames = string(inmatp{:,1}); 
             inmatp(:,1) = [];
             inmatp = inmatp(ccog,:);
             fulls(cidx,:) = table2array(inmat);
             fullp(cidx,:) = table2array(inmatp);
             full_lab{cidx} = append(cblock,'_',cctrl,'_',ccog);
             cidx = cidx + 1;
        end   
    end   
end

%Threshold original for sanity check.
nofdr_fulls = fulls;
nofdr_fulls(fullp >= 0.05) = NaN;

%Save each matrix.
outlist = {fulls,fullp,nofdr_fulls};
outlabs = {'fulls','fullp','fulls_nofdr'};
nout = size(outlist,2);
for oidx=1:nout
    cout = outlist{oidx};
    clab = outlabs{oidx};
    ctab = array2table(cout,'RowNames',string(full_lab),'VariableNames',statelabs);
    outfile = append(outpath,clab,'.csv');
    writetable(ctab,outfile,'WriteRowNames',true);
end

%Vectorize.
[nx,ny] = size(fullp);
pvec = reshape(fullp,[(nx*ny) 1]);

%FDR correct for every test in the experiment.
fdrvec = mafdr(pvec,'BHFDR',1);
fdr_fullp = reshape(fdrvec,[nx ny]);
fdr_fulls = fulls;
fdr_fulls(fdr_fullp >= 0.05) = NaN;

%Save each FDR-corrected matrix.
outlist = {fdr_fulls,fdr_fullp};
outlabs = {'fulls_fdr','fullp_fdr'};
nout = size(outlist,2);
for oidx=1:nout
    cout = outlist{oidx};
    clab = outlabs{oidx};
    ctab = array2table(cout,'RowNames',string(full_lab),'VariableNames',statelabs);
    outfile = append(outpath,clab,'.csv');
    writetable(ctab,outfile,'WriteRowNames',true);
end
disp('Saved.')
end
