function getVoxPref_symmetry_new(isub, visualRegions, sym_type)
% getVoxPref_symmetry_new(isub, visualRegions, sym_type)
%
% Adapted from: Roth, ZN, Kay, K, and Merriam, EP (2022).
% purpose: Extract unique variance gains (Delta R2) and average cross-validation 
%          splits for the hierarchical symmetry models.
%
% Input:
%   isub: Subject index (1-8)
%   visualRegions: Array of regions to process (e.g., 1:8)
%   sym_type: 'mir', 'par', or 'tap'

% --- Setup Paths and Constants --- %
resultsDir = '/bwdata/NSDData/Seohee/Symmetry/regress_results/';
saveDir = '/bwdata/NSDData/Seohee/Symmetry/regress_results/';

sym_methods = {'contour', 'medialAxis', 'area'};
nsplits = 2; % Matches regressPrfSplit

% Setup ROI metadata for saving
visualRoisfolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC', 'LO1', 'LO2'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC','LO'};

% Reconstruct the exact ROI mask used during regression
visualRoisFile = fullfile(visualRoisfolder,'roi/prf-visualrois.nii.gz');
visRoiData = niftiread(visualRoisFile);

placesRoisFile = fullfile(visualRoisfolder,'roi/floc-places.nii.gz'); 
if exist(placesRoisFile, 'file')
    placeRoiData = niftiread(placesRoisFile);
    visRoiData(placeRoiData == 1) = 8;
    visRoiData(placeRoiData == 2) = 9;
    visRoiData(placeRoiData == 3) = 10;
end

kastnerFile = fullfile(visualRoisfolder,'roi/Kastner2015.nii.gz');
if exist(kastnerFile, 'file')
    kastnerData = niftiread(kastnerFile);
    visRoiData(visRoiData == 0 & kastnerData == 15) = 11; % LO1
    visRoiData(visRoiData == 0 & kastnerData == 14) = 12; % LO2
end

% Initialize containers for the final output
allMethodResults = struct();

for m = 1:length(sym_methods)
    curMethod = sym_methods{m};
    disp(['Processing Method: ', curMethod]);
    
    for iregion = visualRegions
        % Load the regression results
        fname = fullfile(resultsDir, sprintf('regressPrfSplit_%s_%s_v%d_sub%d.mat', ...
            sym_type, curMethod, iregion, isub));
        
        if ~exist(fname, 'file')
            warning('File not found: %s', fname);
            continue;
        end
        
        dat = load(fname);
        nsd = dat.nsd;
        rois = dat.rois;
        
        % --- 1. CONCATENATE ROIs (e.g., V1v + V1d -> V1) --- %
        % Regions 1, 2, 3, and 8 have multiple ROIs (ventral/dorsal or LO1/2)
        if length(rois) > 1
            combined = struct();
            combined.M0_r2 = [];
            combined.M1_r2 = [];
            combined.M2_r2 = [];
            combined.F_M2_vs_M1 = [];
            combined.p_M2_vs_M1 = [];
            combined.F_M1_vs_M0 = []; 
            combined.p_M1_vs_M0 = [];
            combined.roiInd = [];
            
            for iroi = 1:length(rois)
                combined.M0_r2 = cat(2, combined.M0_r2, nsd.M0_r2{iroi});
                combined.M1_r2 = cat(2, combined.M1_r2, nsd.M1_r2{iroi});
                combined.M2_r2 = cat(2, combined.M2_r2, nsd.M2_r2{iroi});
                combined.F_M1_vs_M0 = cat(1, combined.F_M1_vs_M0, nsd.F_M1_vs_M0{iroi});
                combined.p_M1_vs_M0 = cat(1, combined.p_M1_vs_M0, nsd.p_M1_vs_M0{iroi});
                combined.F_M2_vs_M1 = cat(1, combined.F_M2_vs_M1, nsd.F_M2_vs_M1{iroi});
                combined.p_M2_vs_M1 = cat(1, combined.p_M2_vs_M1, nsd.p_M2_vs_M1{iroi});
                combined.roiInd = cat(1, combined.roiInd, nsd.roiInd{iroi});
            end
            nsd = combined;
        else
            nsd.M0_r2 = nsd.M0_r2{1};
            nsd.M1_r2 = nsd.M1_r2{1};
            nsd.M2_r2 = nsd.M2_r2{1};
            nsd.F_M1_vs_M0 = nsd.F_M1_vs_M0{1};
            nsd.p_M1_vs_M0 = nsd.p_M1_vs_M0{1};
            nsd.F_M2_vs_M1 = nsd.F_M2_vs_M1{1};
            nsd.p_M2_vs_M1 = nsd.p_M2_vs_M1{1};
            nsd.roiInd = nsd.roiInd{1};
        end
        
        % --- 2. CALCULATE UNIQUE VARIANCE (DELTA R2) --- %
        % Formula: Gain = Model(N) - Model(N-1)
        % Note: These are [nsplits x nvox]
        deltaR2_Skel = nsd.M1_r2 - nsd.M0_r2;
        deltaR2_Sym  = nsd.M2_r2 - nsd.M1_r2;
        totalR2_Full = nsd.M2_r2;
        
        % --- 3. AVERAGE ACROSS SPLITS --- %
        % Store the mean in nsplits + 1 row
        avgM0R2      = mean(nsd.M0_r2, 1);
        avgDeltaSkel = mean(deltaR2_Skel, 1);
        avgDeltaSym  = mean(deltaR2_Sym, 1);
        avgTotalR2   = mean(totalR2_Full, 1);
        
        % Store averaged results for this region and method
        allMethodResults.(curMethod).M0R2{iregion}      = [nsd.M0_r2; avgM0R2];
        allMethodResults.(curMethod).deltaSkel{iregion} = [deltaR2_Skel; avgDeltaSkel];
        allMethodResults.(curMethod).deltaSym{iregion}  = [deltaR2_Sym; avgDeltaSym];
        allMethodResults.(curMethod).totalR2{iregion}   = [totalR2_Full; avgTotalR2];
        allMethodResults.(curMethod).FstatSkel{iregion} = nsd.F_M1_vs_M0;
        allMethodResults.(curMethod).pvalSkel{iregion}  = nsd.p_M1_vs_M0;
        allMethodResults.(curMethod).FstatSym{iregion}  = nsd.F_M2_vs_M1;
        allMethodResults.(curMethod).pvalSym{iregion}   = nsd.p_M2_vs_M1;
        allMethodResults.(curMethod).roiInd{iregion}    = nsd.roiInd;
    end
end

% --- 4. SAVE SUMMARY DATA --- %
% This file will be used for cortical surface projection and ROI plotting
saveFile = fullfile(saveDir, sprintf('voxModelPref_%s_sub%d.mat', sym_type, isub));
save(saveFile, 'allMethodResults', 'combinedRoiNames', 'roiNames', 'sym_methods', 'visRoiData', 'nsplits');
disp(['Summary saved to: ', saveFile]);

end