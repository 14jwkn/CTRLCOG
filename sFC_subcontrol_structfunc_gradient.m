%{
Across participants, read in their sFC matrices, Fisher Z-transform,
average, reverse transform, and calculate gradients. Save the gradients,
and check the variance and first three gradients.
Output:
sFC_gradients.csv Scores for each gradient derived from the average sFC.
sFC_gradients_variance.csv Variance for each gradient derived from the average sFC.
sFC_gradients_variance.jpg Variance for each gradient derived from the average sFC plotted.
sFC_gradients_123.png Scores for the first three gradients derived from the average sFC plotted.
%}

%Define command line arguments.
function [] = sFC_subcontrol_structfunc_gradient(threstype,thresval)
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

%Set parameters.
pkeep = 1 - (str2num(thresval)/100);

%Read in subjects.
subfile = 'dr_full_intersect.txt';       
subjects = string(textread(subfile,'%s','delimiter',','));
nsub = size(subjects,1);

%Read in sFC.
infile = append(sFCpath,'sub_sFC.h5');
inkey = append('/sFC');
sFCmat = h5read(infile,inkey);

%Fisher transform each row.
sFCfish = arrayfun(@(idx) atanh(sFCmat(idx,:)),1:size(sFCmat,1),'UniformOutput', false);
sFCfish = cat(1,sFCfish{:});

%Average.
sFCavg = mean(sFCfish,1);

%Generate square matrix with zeros in diagonal.
nroi = 360;
sq_mat = zeros(nroi,nroi);
sq_mat(triu(true(nroi),1)) = sFCavg;
sq_mat = sq_mat + sq_mat.';

%No threshold.
if strcmp(threstype,'none')

    %Set again.
    A = sq_mat;
end    

%Fit with the same parameters as in Parkes et al. (2021).
gm = GradientMaps('n_components',5,...
                  'approach','dm',...
                  'kernel','normalized_angle',...
                  'random_state',0);
gm = gm.fit(A); 
cgrad = gm.gradients{1};
clam = gm.lambda{1};

%Save matrices.
ctab = cgrad;
outtab = array2table(ctab);
outfile = append(outpath,'sFC_gradients.csv');
writetable(outtab,outfile,'WriteRowNames',0,'WriteVariableNames',0)
ctab = clam;
outtab = array2table(ctab);
outfile = append(outpath,'sFC_gradients_variance.csv');
writetable(outtab,outfile,'WriteRowNames',0,'WriteVariableNames',0)

%Read in flippers to flip the gradient dimension in a desired direction if it exists.
infile = append(outpath,'sFC_gradients_flip.csv');
if isfile(infile)
    yflip = readtable(infile,'readVariableNames',0);
    nflip = size(yflip,2);
    for flidx=1:nflip
        cflip = yflip{1,flidx};
        cflip = cflip{1};
        if strcmp(cflip,'T')
            cgrad(:,flidx) = -cgrad(:,flidx);
        end
    end
end

%Save plots.
figure;
plot(clam)
outfile = append(outpath,'sFC_gradients_variance.jpg');
saveas(gcf,outfile)
gradient_in_euclidean(cgrad(:,1:2))
xlim([min(cgrad(:,1))*1.1 max(cgrad(:,1))*1.1])
ylim([min(cgrad(:,2))*1.1 max(cgrad(:,2))*1.1])
outfile = append(outpath,'sFC_gradients_123.png');
print(gcf,outfile,'-dpng','-r720')
end
