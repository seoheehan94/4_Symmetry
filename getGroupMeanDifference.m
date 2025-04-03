clear all;
method = 'contour'; % Choose from: 'contour', 'medialAxis', 'area'
symTypes = {'par', 'mir', 'tap'};

filedir_base = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/';

addpath(genpath('/usr/local/freesurfer/7.4.1/matlab'));
[~, M, mr] = load_mgh([filedir_base, 'surfaceData_par/sub1/contourBrain_sub1_lh_fsaverage.mgh']);

% Cutoff thresholds for R² filtering
cutoff_lh = 0;
cutoff_rh = 0;

% Initialize storage
vol_lh = struct(); vol_rh = struct();
R2_lh = struct(); R2_rh = struct();
weighted_diff_lh = struct(); weighted_diff_rh = struct();

% Load data for all symmetry types
for itype = 1:length(symTypes)
    type = symTypes{itype};
    filedir = [filedir_base, 'surfaceData_', type, '/'];

    for isub = 1:8
        % Load left and right hemisphere volume data
        vol_lh.(type)(:, isub) = load_mgh([filedir, 'sub', num2str(isub), '/', method, 'Brain_sub', num2str(isub), '_lh_fsaverage.mgh']);
        vol_rh.(type)(:, isub) = load_mgh([filedir, 'sub', num2str(isub), '/', method, 'Brain_sub', num2str(isub), '_rh_fsaverage.mgh']);

        % Load R² values for weighting
        R2_lh.(type)(:, isub) = load_mgh([filedir, 'sub', num2str(isub), '/', method, 'BrainR2_sub', num2str(isub), '_lh_fsaverage.mgh']);
        R2_rh.(type)(:, isub) = load_mgh([filedir, 'sub', num2str(isub), '/', method, 'BrainR2_sub', num2str(isub), '_rh_fsaverage.mgh']);
    end
end

% Compute mean R² across subjects for filtering
meanR2_lh = struct();
meanR2_rh = struct();
for itype = 1:length(symTypes)
    type = symTypes{itype};
    meanR2_lh.(type) = mean(R2_lh.(type), 2, "omitnan");
    meanR2_rh.(type) = mean(R2_rh.(type), 2, "omitnan");
end

% Compute weighted pairwise differences between symmetry types
pairNames = {'par_mir', 'par_tap', 'mir_tap'};
pairs = {'par', 'mir'; 'par', 'tap'; 'mir', 'tap'};

for ipair = 1:length(pairs)
    type1 = pairs{ipair, 1};
    type2 = pairs{ipair, 2};

    % Compute weighted difference voxel-wise for all subjects
    weighted_diff_lh.(pairNames{ipair}) = compute_weighted_difference(vol_lh.(type1), vol_lh.(type2), R2_lh.(type1), R2_lh.(type2));
    weighted_diff_rh.(pairNames{ipair}) = compute_weighted_difference(vol_rh.(type1), vol_rh.(type2), R2_rh.(type1), R2_rh.(type2));

    % Compute final mean across subjects
    mean_weighted_diff_lh = nanmean(weighted_diff_lh.(pairNames{ipair}), 2);
    mean_weighted_diff_rh = nanmean(weighted_diff_rh.(pairNames{ipair}), 2);

    % Apply R² filtering: Set to NaN if mean R² ≤ 0
    max_threshold = 5; % Adjust based on distribution
    min_threshold = -5;

    R2cutoff_diff_lh = mean_weighted_diff_lh;
    R2cutoff_diff_rh = mean_weighted_diff_rh;
    R2cutoff_diff_lh(mean_weighted_diff_lh > max_threshold) = NaN;
    R2cutoff_diff_lh(mean_weighted_diff_lh < min_threshold) = NaN;
    R2cutoff_diff_rh(mean_weighted_diff_rh > max_threshold) = NaN;
    R2cutoff_diff_rh(mean_weighted_diff_rh < min_threshold) = NaN;

    R2cutoff_diff_lh(meanR2_lh.(type1) <= cutoff_lh | meanR2_lh.(type2) <= cutoff_lh) = NaN;
    R2cutoff_diff_rh(meanR2_rh.(type1) <= cutoff_rh | meanR2_rh.(type2) <= cutoff_rh) = NaN;


    % Save full & filtered results
    filedir_save = [filedir_base, 'surfaceData_diff/'];
    save_mgh(mean_weighted_diff_lh, [filedir_save, method, '_lh_', pairNames{ipair}, '_weighted_diff_full.mgh'], M, mr);
    save_mgh(mean_weighted_diff_rh, [filedir_save, method, '_rh_', pairNames{ipair}, '_weighted_diff_full.mgh'], M, mr);

    save_mgh(R2cutoff_diff_lh, [filedir_save, method, '_lh_', pairNames{ipair}, '_weighted_diff_R2cut.mgh'], M, mr);
    save_mgh(R2cutoff_diff_rh, [filedir_save, method, '_rh_', pairNames{ipair}, '_weighted_diff_R2cut.mgh'], M, mr);
end

%% **Function for Weighted Difference Calculation**
function weighted_diff = compute_weighted_difference(vol1, vol2, R2_1, R2_2)
    % Ensure R² values are non-negative
    R2_1(R2_1 < 0) = 0;
    R2_2(R2_2 < 0) = 0;

    % Compute weighted difference
    weighted_diff = nan(size(vol1));
    for ivox = 1:size(vol1, 1)
        validMask = ~isnan(vol1(ivox, :)) & ~isnan(vol2(ivox, :)) & ~isnan(R2_1(ivox, :)) & ~isnan(R2_2(ivox, :));
        
        if any(validMask)
            % Weighted difference formula
            numerator = R2_1(ivox, validMask) .* vol1(ivox, validMask) - R2_2(ivox, validMask) .* vol2(ivox, validMask);
            denominator = R2_1(ivox, validMask) + R2_2(ivox, validMask);
            denominator(denominator == 0) = NaN; % Avoid division by zero

            weighted_diff(ivox, validMask) = numerator ./ denominator;
        end
    end
end
