
% regressPrfSplit_symmetry_.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: regressPrfSplit(1,1)
%   by: zvi roth
%   date: 7/29/2022
%   purpose: Perform linear regression on the response amplitudes for each voxel with filter output values as predictors
%   uses files created by: prfSampleModel.m, prfSampleModel_synth.m
%   creates files used by: getVoxPref.m

% sym_type: 'mir', 'tap', 'par'
% sym_method: 'contour', 'medialAxis', 'area'

function regressPrfSplit_symmetry_new(isub,visualRegions,sym_type)
% This script processes all three spatial methods: 'contour', 'medialAxis', 'area'

clearvars -except isub visualRegions sym_type
addpath('/home/hanseohe/Documents/GitHub/2_Orientation_Tuning/EXP2/model_computation')

tic

%%%%%%%%%%%%%%%%%%%%%%%%
nsessionsSub = [40 40 32 30 40 32 40 30];
nsessions=nsessionsSub(isub);
nsplits=2;

% Define all spatial methods to loop through
sym_methods = {'contour', 'medialAxis', 'area'};

% Define directories
boxfolder_contour = '/bwdata/NSDData/Seohee/Symmetry/prfsample_ma/';
boxfolder_skeleton = '/bwdata/NSDData/Seohee/Symmetry/prfsample_ma/';
boxfolder_symmetry = ['/bwdata/NSDData/Seohee/Symmetry/prfsample_',sym_type,'/'];
betasfolder = ['/bwdata/NSDData/nsddata_betas/ppdata/subj0' num2str(isub) '/func1pt8mm/betas_fithrf_GLMdenoise_RR/'];
nsdfolder = '/bwdata/NSDData/nsddata/experiments/nsd/';
visualRoisfolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];

% --- ROI Setup --- %
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC', 'LO1', 'LO2'};
visualRoisFile = fullfile(visualRoisfolder,'roi/prf-visualrois.nii.gz');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);

placesRoisFile = fullfile(visualRoisfolder,'roi/floc-places.nii.gz'); %OPA, PPA, RSC 
placeRoiData = niftiread(placesRoisFile);

kastnerFile = fullfile(visualRoisfolder,'roi/Kastner2015.nii.gz');
kastnerData = niftiread(kastnerFile);

visRoiData(placeRoiData == 1) = 8;
visRoiData(placeRoiData == 2) = 9;
visRoiData(placeRoiData == 3) = 10;
visRoiData(visRoiData == 0 & kastnerData == 15) = 11; % Label LO1 as 11
visRoiData(visRoiData == 0 & kastnerData == 14) = 12; % Label LO2 as 12
visRoiData = visRoiData(:);

for visualRegion=visualRegions
    disp(['Processing Visual Region: ', num2str(visualRegion)]);

    %% 1. LOAD DATA
    % Load Contours
    load(fullfile(boxfolder_contour,['prfSampleStim_ma_contour_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),...
        'prfSampleLevOri','rois','allImgs');
    prf_contour = prfSampleLevOri;

    % Load Medial Axis (Skeleton)
    load(fullfile(boxfolder_skeleton,['prfSampleStim_ma_medialAxis_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']), 'prfSampleLevOri');
    prf_skeleton = prfSampleLevOri;
    
    % Load Betas
    for roinum=1:length(rois); iroi = rois(roinum); roiBetas{roinum}=[]; end
    for isession=1:nsessions
        betasfilename = fullfile(betasfolder,['betas_session' num2str(isession,'%02.f') '.nii.gz']);
        betas = niftiread(betasfilename);
        betas = cast(betas,'double');
        betas = betas/300;
        betas=reshape(betas,[],size(betas,4));
        for roinum=1:length(rois)
            iroi = rois(roinum);
            roiBetas{roinum} = [roiBetas{roinum} betas(visRoiData==iroi,:)];
            roiInd{roinum} = find(visRoiData==iroi);
        end
    end
    
    % Setup Trials and Splits
    nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
    nsdDesign = load(nsdDesignFilename);
    nsdDesign.masterordering;%for each of 30000 trials, what is corresponding image (out of 10000 images)
    subDesign = nsdDesign.subjectim(isub,nsdDesign.masterordering);%for each of 30000 trials, what is corresponding image (out of 73000 images)
    [imgTrials, imgNum] = ismember(subDesign, allImgs);%logical array
    
    %if less than 40 sessions, only use image trials that were actually presented
    imgTrials = imgTrials(1:size(roiBetas{roinum},2));
    imgNum = imgNum(1:size(roiBetas{roinum},2));
    splitImgTrials = repmat(imgTrials,2,1);
    midImg = ceil(median(imgNum));
    splitImgTrials(1,imgNum<midImg) = zeros;
    splitImgTrials(2,imgNum>=midImg) = zeros;

    % Identify all valid trials for the full-dataset F-test
    allValidTrials = (imgTrials > 0);
    n_total_trials = sum(allValidTrials);
  
    
%% 2. INNER LOOP: PROCESS EACH SYMMETRY METHOD
    for m = 1:length(sym_methods)
        sym_method = sym_methods{m};
        disp(['--- Running regressions for method: ', sym_method, ' ---']);
        
        % Load specific Symmetry predictor
        load(fullfile(boxfolder_symmetry,['prfSampleStim_', sym_type, '_', sym_method,'_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']), 'prfSampleLevOri');
        prf_symmetry = prfSampleLevOri;
        
        nsd = struct(); % Initialize output struct for this method
        nsd.roiInd = roiInd;
        
        for roinum=1:length(rois)
            nvox = size(roiBetas{roinum},1);
            
            % Trackers for Split-Half CV R2
            r2_split_M0 = zeros(nsplits, nvox);
            r2_split_M1 = zeros(nsplits, nvox);
            r2_split_M2 = zeros(nsplits, nvox);
            
            % Trackers for Partial F-Test Statistics (Full Dataset)
            F_M1_vs_M0 = zeros(nvox, 1); p_M1_vs_M0 = zeros(nvox, 1);
            F_M2_vs_M1 = zeros(nvox, 1); p_M2_vs_M1 = zeros(nvox, 1);
            
            % Predictor counts (including constant)
            p0 = 2; % Contour + const
            p1 = 3; % Contour + Skeleton + const
            p2 = 4; % Contour + Skeleton + Symmetry + const

            %% A. Cross-Validated R2 (Split-Half)
            coef_M0 = zeros(nsplits, nvox, p0);
            coef_M1 = zeros(nsplits, nvox, p1);
            coef_M2 = zeros(nsplits, nvox, p2);
            
            % 1. Train on Splits
            for isplit=1:nsplits
                imgT = splitImgTrials(isplit,:);
                valid_imgs = imgNum(imgT>0);
                for ivox=1:nvox
                    voxBetas = roiBetas{roinum}(ivox, imgT>0)';
                    
                    X0 = [squeeze(prf_contour{roinum}(valid_imgs, ivox, 1, 1)), ones(length(valid_imgs), 1)];
                    X1 = [X0(:,1), squeeze(prf_skeleton{roinum}(valid_imgs, ivox, 1, 1)), X0(:,2)];
                    X2 = [X1(:,1:2), squeeze(prf_symmetry{roinum}(valid_imgs, ivox, 1, 1)), X1(:,3)];
                    
                    coef_M0(isplit, ivox, :) = X0 \ voxBetas;
                    coef_M1(isplit, ivox, :) = X1 \ voxBetas;
                    if rank(X2) == size(X2, 2), coef_M2(isplit, ivox, :) = X2 \ voxBetas; else, coef_M2(isplit, ivox, :) = NaN; end
                end
            end
            
            % 2. Test on Opposite Split
            for isplit=1:nsplits
                imgT = splitImgTrials(isplit,:);
                valid_imgs = imgNum(imgT>0);
                other_split = nsplits - isplit + 1; 
                
                for ivox=1:nvox
                    voxBetas = roiBetas{roinum}(ivox, imgT>0)';
                    
                    X0 = [squeeze(prf_contour{roinum}(valid_imgs, ivox, 1, 1)), ones(length(valid_imgs), 1)];
                    X1 = [X0(:,1), squeeze(prf_skeleton{roinum}(valid_imgs, ivox, 1, 1)), X0(:,2)];
                    X2 = [X1(:,1:2), squeeze(prf_symmetry{roinum}(valid_imgs, ivox, 1, 1)), X1(:,3)];
                    
                    pred_M0 = X0 * squeeze(coef_M0(other_split, ivox, :));
                    pred_M1 = X1 * squeeze(coef_M1(other_split, ivox, :));
                    pred_M2 = X2 * squeeze(coef_M2(other_split, ivox, :));
                    
                    r2_split_M0(isplit, ivox) = rsquared(voxBetas - pred_M0, voxBetas);
                    r2_split_M1(isplit, ivox) = rsquared(voxBetas - pred_M1, voxBetas);
                    r2_split_M2(isplit, ivox) = rsquared(voxBetas - pred_M2, voxBetas);
                end
            end
            
            %% B. Partial F-Tests (Full Dataset for maximum statistical power)
            valid_imgs_all = imgNum(allValidTrials);
            
            for ivox=1:nvox
                voxBetas_all = roiBetas{roinum}(ivox, allValidTrials)';
                
                X0 = [squeeze(prf_contour{roinum}(valid_imgs_all, ivox, 1, 1)), ones(n_total_trials, 1)];
                X1 = [X0(:,1), squeeze(prf_skeleton{roinum}(valid_imgs_all, ivox, 1, 1)), X0(:,2)];
                X2 = [X1(:,1:2), squeeze(prf_symmetry{roinum}(valid_imgs_all, ivox, 1, 1)), X1(:,3)];
                
                % Fit full models
                b0 = X0 \ voxBetas_all; 
                b1 = X1 \ voxBetas_all;
                
                % Calculate Residual Sum of Squares (RSS)
                RSS0 = sum((voxBetas_all - (X0 * b0)).^2);
                RSS1 = sum((voxBetas_all - (X1 * b1)).^2);
                
                % F-test: Model 1 vs Model 0 (Does Skeleton add value?)
                F_M1_vs_M0(ivox) = ((RSS0 - RSS1) / (p1 - p0)) / (RSS1 / (n_total_trials - p1));
                p_M1_vs_M0(ivox) = 1 - fcdf(F_M1_vs_M0(ivox), (p1 - p0), (n_total_trials - p1));
                
                % F-test: Model 2 vs Model 1 (Does Symmetry add value?)
                if rank(X2) == size(X2, 2)
                    b2 = X2 \ voxBetas_all;
                    RSS2 = sum((voxBetas_all - (X2 * b2)).^2);
                    
                    F_M2_vs_M1(ivox) = ((RSS1 - RSS2) / (p2 - p1)) / (RSS2 / (n_total_trials - p2));
                    p_M2_vs_M1(ivox) = 1 - fcdf(F_M2_vs_M1(ivox), (p2 - p1), (n_total_trials - p2));
                else
                    F_M2_vs_M1(ivox) = NaN; p_M2_vs_M1(ivox) = NaN;
                end
            end
            
            % Save to struct
            nsd.M0_r2{roinum} = r2_split_M0;
            nsd.M1_r2{roinum} = r2_split_M1;
            nsd.M2_r2{roinum} = r2_split_M2;
            nsd.F_M1_vs_M0{roinum} = F_M1_vs_M0; 
            nsd.p_M1_vs_M0{roinum} = p_M1_vs_M0;
            nsd.F_M2_vs_M1{roinum} = F_M2_vs_M1; 
            nsd.p_M2_vs_M1{roinum} = p_M2_vs_M1;
        end
        
        %% SAVE RESULTS FOR THIS METHOD
        savefolder = '/bwdata/NSDData/Seohee/Symmetry/regress_results/';
        if ~exist(savefolder, 'dir'), mkdir(savefolder); end
        
        save(fullfile(savefolder,['regressPrfSplit_', sym_type, '_', sym_method, '_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']), 'nsd', 'rois');
    end
    toc
end
end