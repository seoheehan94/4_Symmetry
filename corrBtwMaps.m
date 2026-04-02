close all;
clear all;

% --- 1. Define Directories ---
parfolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/parfilter/';%to save model outputs
parList = dir(parfolder);
parList = parList(~ismember({parList(:).name},{'.','..'}));

mirfolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/mirfilter/';%to save model outputs
mirList = dir(mirfolder);
mirList = mirList(~ismember({mirList(:).name},{'.','..'}));

tapfolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/mirfilter/';%to save model outputs
tapList = dir(tapfolder);
tapList = tapList(~ismember({tapList(:).name},{'.','..'}));

savefolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/corrBtwMaps/';%to save model outputs
filePath = fullfile(savefolder, 'results_corrBtwMaps_nonzeros.mat');
if isfile(filePath)
    load(filePath, 'allResultsTable');
end
%% --- 2. Process Correlations ---

iimg=0;
for imgNum=1:length(parList)
    iimg = iimg+1
    parname = parList(imgNum).name;
    mirname = mirList(imgNum).name;
    tapname = tapList(imgNum).name;


    parModel = load(fullfile(parfolder, parname), 'model');
    mirModel = load(fullfile(mirfolder, mirname), 'model');
    tapModel = load(fullfile(tapfolder, tapname), 'model');

    curpar=parModel.model.contour;
    curmir=mirModel.model.contour;
    curtap=tapModel.model.contour;

    % --- Streamlined Masking Logic ---
    % Find indices where AT LEAST ONE of the maps is non-zero.
    % This excludes background (mutual zeros) but keeps a 0 if one map 
    % has a signal and the others do not.
    valid_idx = (curpar ~= 0) | (curmir ~= 0) | (curtap ~= 0);
    
    curpar_valid = curpar(valid_idx);
    curmir_valid = curmir(valid_idx);
    curtap_valid = curtap(valid_idx);


    % --- Calculate 3 Pairwise Correlations ---
    % 1. Parallelism vs Mirror
    [R_pm, P_pm] = corrcoef(curpar_valid, curmir_valid);
    
    % 2. Parallelism vs Taper
    [R_pt, P_pt] = corrcoef(curpar_valid, curtap_valid);
    
    % 3. Mirror vs Taper
    [R_mt, P_mt] = corrcoef(curmir_valid, curtap_valid);

    % --- Store Results ---
    % Extracting the off-diagonal elements (1,2) from the matrices
    resultsTable = table(imgNum, R_pm(1,2), P_pm(1,2), R_pt(1,2), P_pt(1,2), R_mt(1,2), P_mt(1,2), ...
        'VariableNames', {'imgnum', 'R_par_mir', 'P_par_mir', 'R_par_tap', 'P_par_tap', 'R_mir_tap', 'P_mir_tap'});
    
    if exist('allResultsTable', 'var')
        allResultsTable = [allResultsTable; resultsTable];
    else
        allResultsTable = resultsTable;
    end
    
    save(filePath, 'allResultsTable');

  
end

writetable(allResultsTable, 'results_corrBtwMaps.csv');