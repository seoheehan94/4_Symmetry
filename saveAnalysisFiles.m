clear all;

% --- Parameters --- %
sym_type = 'tap'; % sym_type = {'mir', 'tap', 'par'};
sym_methods = {'contour', 'medialAxis', 'area'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC','LO'};

prffolder = '/bwdata/NSDData/Seohee/Symmetry/regress_results/';
save_brainfolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/brainVolume/';
save_analysisfolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/analysis/';
if ~exist(save_brainfolder, 'dir'), mkdir(save_brainfolder); end
if ~exist(save_analysisfolder, 'dir'), mkdir(save_analysisfolder); end

%% 1. CREATE CSV FOR R ANALYSIS %%
disp('Compiling data for CSV...');

subj_all = [];
roi_all  = {};
method_all = {};
M0R2_all = [];
deltaSkel_all = [];
deltaSym_all  = [];
totalR2_all = [];
FstatSkel_all = [];
FstatSym_all = [];  
pvalSkel_all = [];
pvalSym_all = [];

for isub = 1:8
    fprintf('Loading Subject %d...\n', isub);
    fname = fullfile(prffolder, sprintf('voxModelPref_%s_sub%d.mat', sym_type, isub));
    if ~exist(fname, 'file'), continue; end
    
    load(fname); % Loads: allMethodResults, nsplits, visRoiData
    avgIdx = nsplits + 1; % The row containing the average across splits
    
    for m = 1:length(sym_methods)
        curMethod = sym_methods{m};
        
        for r = 1:length(combinedRoiNames)
            roiName = combinedRoiNames{r};
            
            % Extract data for this ROI
            % (Check if region exists, occasionally some subjects miss a top-level ROI)
            if length(allMethodResults.(curMethod).roiInd) < r || isempty(allMethodResults.(curMethod).roiInd{r})
                continue; 
            end
            
            % Get the averaged values (last row) and stats
            tmp_M0R2      = allMethodResults.(curMethod).M0R2{r}(avgIdx, :);
            tmp_deltaSkel = allMethodResults.(curMethod).deltaSkel{r}(avgIdx, :);
            tmp_deltaSym  = allMethodResults.(curMethod).deltaSym{r}(avgIdx, :);
            tmp_totalR2   = allMethodResults.(curMethod).totalR2{r}(avgIdx, :);
            
            tmp_FstatSkel = allMethodResults.(curMethod).FstatSkel{r}; 
            tmp_FstatSym  = allMethodResults.(curMethod).FstatSym{r};  
            tmp_pvalSkel  = allMethodResults.(curMethod).pvalSkel{r};
            tmp_pvalSym   = allMethodResults.(curMethod).pvalSym{r};

            tmp_M0R2      = tmp_M0R2(:);
            tmp_deltaSkel = tmp_deltaSkel(:);
            tmp_deltaSym  = tmp_deltaSym(:);
            tmp_totalR2   = tmp_totalR2(:);
            tmp_FstatSkel = tmp_FstatSkel(:);
            tmp_FstatSym  = tmp_FstatSym(:);
            tmp_pvalSkel  = tmp_pvalSkel(:);
            tmp_pvalSym   = tmp_pvalSym(:);
            nvox = numel(tmp_M0R2);
            
            % Append to long-format arrays
            subj_all      = [subj_all; repmat(isub, nvox, 1)];
            roi_all       = [roi_all; repmat({roiName}, nvox, 1)];
            method_all    = [method_all; repmat({curMethod}, nvox, 1)];
            
            M0R2_all      = [M0R2_all; tmp_M0R2];
            deltaSkel_all = [deltaSkel_all; tmp_deltaSkel];
            deltaSym_all  = [deltaSym_all; tmp_deltaSym];
            totalR2_all   = [totalR2_all; tmp_totalR2];
            
            FstatSkel_all = [FstatSkel_all; tmp_FstatSkel]; 
            FstatSym_all  = [FstatSym_all; tmp_FstatSym];   
            pvalSkel_all  = [pvalSkel_all; tmp_pvalSkel];
            pvalSym_all   = [pvalSym_all; tmp_pvalSym];
        end
    end
end

% Build and save table
T_all = table(subj_all, roi_all, method_all, M0R2_all, deltaSkel_all, deltaSym_all, totalR2_all, ...
              FstatSkel_all, FstatSym_all, pvalSkel_all, pvalSym_all, ...
             'VariableNames', {'subj', 'ROI', 'Method', 'M0_R2', 'DeltaR2_Skel', 'DeltaR2_Sym', ...
                               'Total_R2', 'F_Skel', 'F_Sym', 'p_val_Skel', 'p_val_Sym'});
csvFile = fullfile(save_analysisfolder, sprintf('allROI_%s_Summary.csv', sym_type));
writetable(T_all, csvFile);
disp(['CSV saved to: ', csvFile]);


%% 2. CREATE NIFTI BRAIN VOLUMES %%
disp('Generating NIfTI volumes...');
for isub = 1:8
    fprintf('isub:%d. ...\n', isub);
    fname = fullfile(prffolder, sprintf('voxModelPref_%s_sub%d.mat', sym_type, isub));
    if ~exist(fname, 'file'), continue; end
    
    load(fname); % Loads: allMethodResults, nsplits, visRoiData
    avgIdx = nsplits + 1;
    
    % Reference NIfTI for headers
    refFile = ['/bwdata/NSDData/nsddata_betas/ppdata/subj0' num2str(isub) '/func1pt8mm/betas_fithrf_GLMdenoise_RR/betas_session01.nii.gz'];
    if ~exist(refFile, 'file')
        warning('Reference file missing for Sub %d. Skipping NIfTI creation.', isub);
        continue;
    end
    info_ref = niftiinfo(refFile);
    volSize  = size(visRoiData); 
    
    %% ========================================================
    %% A. PROCESS SKELETON (M1) - Runs ONCE per subject
    %% ========================================================
    % M1 is identical across M2 methods, so we pull from the first method
    baseMethod = sym_methods{1}; 
    
    brain_deltaSkel = NaN(volSize);
    brain_pvalSkel  = NaN(volSize);
    brain_M0R2      = NaN(volSize);
    brain_M1R2      = NaN(volSize);
    
    for r = 1:length(combinedRoiNames)
        if length(allMethodResults.(baseMethod).roiInd) < r || isempty(allMethodResults.(baseMethod).roiInd{r})
            continue;
        end
        curIndices = allMethodResults.(baseMethod).roiInd{r};
        
        cur_deltaSkel = allMethodResults.(baseMethod).deltaSkel{r}(avgIdx, :);
        cur_M0R2      = allMethodResults.(baseMethod).M0R2{r}(avgIdx, :);
        
        brain_deltaSkel(curIndices) = cur_deltaSkel;
        brain_pvalSkel(curIndices)  = allMethodResults.(baseMethod).pvalSkel{r};
        brain_M0R2(curIndices)      = cur_M0R2;
        brain_M1R2(curIndices)      = cur_M0R2 + cur_deltaSkel; 
    end
    
    % --- Skeleton Masks ---
    mask_Skel_Prev = (brain_pvalSkel < 0.05) & (brain_deltaSkel > 0);
    brain_deltaSkel_PrevMasked = brain_deltaSkel;
    brain_deltaSkel_PrevMasked(~mask_Skel_Prev) = NaN;
    
    mask_Skel_Mag = (brain_M1R2 > 0);
    brain_deltaSkel_MagMasked = brain_deltaSkel;
    brain_deltaSkel_MagMasked(~mask_Skel_Mag) = NaN;
    
    % --- Skeleton Outputs (No 'curMethod' in the names!) ---
    skel_outputs = {
        brain_deltaSkel,            'Brain_Skel_Raw_DeltaR2';
        brain_deltaSkel_PrevMasked, 'Brain_Skel_Prevalence_DualCriterion';
        brain_deltaSkel_MagMasked,  'Brain_Skel_Magnitude_PosM1R2';
        brain_M1R2,                 'Brain_Skel_M1R2';
    };
    
    % --- Write Skeleton NIfTIs ---
    for i = 1:size(skel_outputs,1)
        data = skel_outputs{i,1};
        if length(size(data)) > 3, data = squeeze(data); end
        
        outName = fullfile(save_brainfolder, [skel_outputs{i,2}, '_sub', num2str(isub), '.nii']);
        niftiwrite(data, outName);
        
        info_new = niftiinfo(outName);
        info_new.PixelDimensions = [1.8 1.8 1.8];
        info_new.TransformName     = info_ref.TransformName;
        info_new.SpatialDimension  = info_ref.SpatialDimension;
        info_new.Transform         = info_ref.Transform;
        info_new.Qfactor           = info_ref.Qfactor;
        if isfield(info_ref, 'raw')
            info_new.raw.sform_code = info_ref.raw.sform_code;
            info_new.raw.srow_x     = info_ref.raw.srow_x;
            info_new.raw.srow_y     = info_ref.raw.srow_y;
            info_new.raw.srow_z     = info_ref.raw.srow_z;
            info_new.raw.pixdim     = info_ref.raw.pixdim;
        end
        niftiwrite(data, outName, info_new);
    end

    %% ========================================================
    %% B. PROCESS SYMMETRY (M2) - Runs for EACH spatial method
    %% ========================================================
    for m = 1:length(sym_methods)
        curMethod = sym_methods{m};
        
        brain_deltaSym  = NaN(volSize);
        brain_pvalSym   = NaN(volSize);
        brain_totalR2   = NaN(volSize);
        
        for r = 1:length(combinedRoiNames)
            if length(allMethodResults.(curMethod).roiInd) < r || isempty(allMethodResults.(curMethod).roiInd{r})
                continue;
            end
            curIndices = allMethodResults.(curMethod).roiInd{r};
            
            brain_deltaSym(curIndices) = allMethodResults.(curMethod).deltaSym{r}(avgIdx, :);
            brain_pvalSym(curIndices)  = allMethodResults.(curMethod).pvalSym{r};
            brain_totalR2(curIndices)  = allMethodResults.(curMethod).totalR2{r}(avgIdx, :);
        end
        
        % --- Symmetry Masks ---
        mask_Sym_Prev = (brain_pvalSym < 0.05) & (brain_deltaSym > 0);
        brain_deltaSym_PrevMasked = brain_deltaSym;
        brain_deltaSym_PrevMasked(~mask_Sym_Prev) = NaN;
        
        mask_Sym_Mag = (brain_totalR2 > 0);
        brain_deltaSym_MagMasked = brain_deltaSym;
        brain_deltaSym_MagMasked(~mask_Sym_Mag) = NaN;
        
        % --- Symmetry Outputs ---
        sym_outputs = {
            brain_deltaSym,            sprintf('Brain_Sym_%s_%s_Raw_DeltaR2', sym_type, curMethod);
            brain_deltaSym_PrevMasked, sprintf('Brain_Sym_%s_%s_Prevalence_DualCriterion', sym_type, curMethod);
            brain_deltaSym_MagMasked,  sprintf('Brain_Sym_%s_%s_Magnitude_PosTotalR2', sym_type, curMethod);
            brain_totalR2,             sprintf('Brain_Sym_%s_%s_totalR2', sym_type, curMethod);
        };
        
        % --- Write Symmetry NIfTIs ---
        for i = 1:size(sym_outputs,1)
            data = sym_outputs{i,1};
            if length(size(data)) > 3, data = squeeze(data); end
            
            outName = fullfile(save_brainfolder, [sym_outputs{i,2}, '_sub', num2str(isub), '.nii']);
            niftiwrite(data, outName);
            
            info_new = niftiinfo(outName);
            info_new.PixelDimensions = [1.8 1.8 1.8];
            info_new.TransformName     = info_ref.TransformName;
            info_new.SpatialDimension  = info_ref.SpatialDimension;
            info_new.Transform         = info_ref.Transform;
            info_new.Qfactor           = info_ref.Qfactor;
            if isfield(info_ref, 'raw')
                info_new.raw.sform_code = info_ref.raw.sform_code;
                info_new.raw.srow_x     = info_ref.raw.srow_x;
                info_new.raw.srow_y     = info_ref.raw.srow_y;
                info_new.raw.srow_z     = info_ref.raw.srow_z;
                info_new.raw.pixdim     = info_ref.raw.pixdim;
            end
            niftiwrite(data, outName, info_new);
        end
    end
end
disp('All NIfTI files generated.');