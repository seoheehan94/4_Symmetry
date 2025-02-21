
% regressPrfSplit_length.m
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

function regressPrfSplit_symmetry(isub,visualRegions,type,method)
% addpath(genpath('/home/hanseohe/Documents/GitHub/nsdOtopy'));
addpath('/home/hanseohe/Documents/GitHub/2_Orientation_Tuning/EXP2/model_computation')

tic

%%%%%%%%%%%%%%%%%%%%%%%%
nsessionsSub = [40 40 32 30 40 32 40 30];
nsessions=nsessionsSub(isub);
nsplits=2;
bandpass = 1; bandMin = 1; bandMax = 1;

boxfolder = ['/bwdata/NSDData/Seohee/Symmetry/prfsample_',type,'/'];
betasfolder = ['/bwdata/NSDData/nsddata_betas/ppdata/subj0' num2str(isub) '/func1pt8mm/betas_fithrf_GLMdenoise_RR/'];
% stimfilename = fullfile(folder,'nsdsynthetic_colorstimuli_subj01.hdf5');
nsdfolder = '/bwdata/NSDData/nsddata/experiments/nsd/';
visualRoisfolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];

visualRoisFile = fullfile(visualRoisfolder,'roi/prf-visualrois.nii.gz');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
placesRoisFile = fullfile(visualRoisfolder,'roi/floc-places.nii.gz'); %OPA, PPA, RSC 
placeRoiData = niftiread(placesRoisFile);
visRoiData(placeRoiData == 1) = 8;
visRoiData(placeRoiData == 2) = 9;
visRoiData(placeRoiData == 3) = 10;
visRoiData = visRoiData(:);

for visualRegion=visualRegions
    visualRegion
    load(fullfile(boxfolder,['prfSampleStim_',type '_', method,'_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri',...
        'rois','allImgs','numLevels','numFeatures','interpImgSize','backgroundSize','pixPerDeg',...
        'roiPrf');
    %if prf sampling was done with the nonlinear CSS prf, then we want to
    %define the weights for the constrained model as a sum of the
    %orientation model across orientations:
     for roinum=1:length(rois)
         prfSampleLev{roinum} = squeeze(sum(prfSampleLevOri{roinum},4));
     end
    
    if bandpass
        for roinum=1:length(rois)
            prfSampleLevOri{roinum} = prfSampleLevOri{roinum}(:,:,bandMin:bandMax,:);
             prfSampleLev{roinum} = prfSampleLev{roinum}(:,:,bandMin:bandMax);
        end
        numLevels = bandMax-bandMin+1;
    end
    
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
    
    nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
    nsdDesign = load(nsdDesignFilename);
    nsdDesign.masterordering;%for each of 30000 trials, what is corresponding image (out of 10000 images)
    subDesign = nsdDesign.subjectim(isub,nsdDesign.masterordering);%for each of 30000 trials, what is corresponding image (out of 73000 images)
    [imgTrials, imgNum] = ismember(subDesign, allImgs);%logical array
    
    %if less than 40 sessions, only use image trials that were actually presented
    imgTrials = imgTrials(1:size(roiBetas{roinum},2));
    imgNum = imgNum(1:size(roiBetas{roinum},2));
    
    imgTrialSum = cumsum(imgTrials);
    splitImgTrials = repmat(imgTrials,2,1);
    
    midImg = ceil(median(imgNum));
    splitImgTrials(1,imgNum<midImg) = zeros;
    splitImgTrials(2,imgNum>=midImg) = zeros;
    maxNumTrials = max(sum(splitImgTrials,2));
    
    r2 = cell(length(rois),1);
    r2feat = cell(length(rois),1);
    r2featSplit = cell(length(rois),1);
    r2split = cell(length(rois),1);
    voxFeatResidFeatR2 = cell(length(rois),1);
    voxFeatPredFeatR2 = cell(length(rois),1);
    voxResidFeatR2 = cell(length(rois),1);
    voxPredFeatR2 = cell(length(rois),1);
    aicFeatSplit = cell(length(rois),1);
    bicFeatSplit = cell(length(rois),1);
    
    pearsonRfeat = cell(length(rois),1);
    pearsonR = cell(length(rois),1);
    
    
    for roinum=1:length(rois)
        nvox(roinum) = size(roiBetas{roinum},1);
        voxFeatResidual{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials);
        voxResidual{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials);
        voxFeatResidualSplit{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials);
        voxResidualSplit{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials);
        voxFeatCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numFeatures+1);
        voxCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels+1);
        voxPredFeatCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numFeatures+1);
        voxPredCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels+1);
        voxFeatPredFeatCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numFeatures+1);
        voxResidFeatCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numFeatures+1);
        voxFeatResidFeatCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numFeatures+1);
        
        %get model coefficients for each voxel, within each split
        for isplit=1:nsplits
            imgTrials = splitImgTrials(isplit,:);
            numTrials = sum(imgTrials);
            for ivox=1:nvox(roinum)
                voxBetas = roiBetas{roinum}(ivox,imgTrials>0)';
                voxPrfSample = squeeze(prfSampleLev{roinum}(imgNum(imgTrials>0),ivox,:));
                %add constant predictor
                voxPrfSample(:,end+1) = ones;
                voxCoef{roinum}(isplit,ivox,:) = voxPrfSample\voxBetas;
                
                voxPrfFeatSample = squeeze(prfSampleLevOri{roinum}(imgNum(imgTrials>0),ivox,:,:));
                voxPrfFeatSample = reshape(voxPrfFeatSample,[],numLevels*numFeatures);
                
                %add constant predictor
                voxPrfFeatSample(:,end+1) = ones;

                rank_X = rank(voxPrfFeatSample);
                if rank_X < size(voxPrfFeatSample, 2)
                    fprintf('isplit: %d, ivox: %d\n', isplit, ivox);
                    warning('Design matrix is rank-deficient.');
                    voxFeatCoef{roinum}(isplit,ivox,:) = [NaN;NaN];
                else
                    voxFeatCoef{roinum}(isplit,ivox,:) = voxPrfFeatSample\voxBetas;%check vox 144 in first ROI
                end
                

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %regress vignetting predicted timecourse with orientation model
                voxPred = squeeze(voxCoef{roinum}(isplit,ivox,:))'*voxPrfSample';
                voxPredFeatCoef{roinum}(isplit,ivox,:) = regress(voxPred',voxPrfFeatSample);
                
                voxFeatPred = squeeze(voxFeatCoef{roinum}(isplit,ivox,:))'*voxPrfFeatSample';
                voxFeatPredFeatCoef{roinum}(isplit,ivox,:) = regress(voxFeatPred',voxPrfFeatSample);
                
                %compute within-split residuals
                voxFeatResidual{roinum}(isplit,ivox,1:numTrials) = voxBetas' - squeeze(voxFeatCoef{roinum}(isplit,ivox,:))'*voxPrfFeatSample';
                voxResidual{roinum}(isplit,ivox,1:numTrials) = voxBetas' - squeeze(voxCoef{roinum}(isplit,ivox,:))'*voxPrfSample';
                
                %regress residuals of vignetting model with full model
                voxResidFeatCoef{roinum}(isplit,ivox,:) = regress(squeeze(voxResidual{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:)))),voxPrfFeatSample);
                %regress residuals of orientation model with orientation model
                voxFeatResidFeatCoef{roinum}(isplit,ivox,:) = regress(squeeze(voxFeatResidual{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:)))),voxPrfFeatSample);
                
                %r2 within split
                r2{roinum}(isplit,ivox) = rsquared(voxResidual{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                r2feat{roinum}(isplit,ivox) = rsquared(voxFeatResidual{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                %resid r2 within split
                %residuals of full orientation model:
                featResidFeatResidual = voxBetas' - squeeze(voxFeatResidFeatCoef{roinum}(isplit,ivox,:))'*voxPrfFeatSample';
                voxFeatResidFeatR2{roinum}(isplit,ivox) = rsquared(featResidFeatResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                %prediction of full orientation model:
                featPredFeatResidual = voxBetas' - squeeze(voxFeatPredFeatCoef{roinum}(isplit,ivox,:))'*voxPrfFeatSample';
                voxFeatPredFeatR2{roinum}(isplit,ivox) = rsquared(featPredFeatResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                
                %residuals of constrained model:
                % residOriResidual = voxBetas' - squeeze(voxResidFeatCoef{roinum}(isplit,ivox,:))'*voxPrfFeatSample';
                % voxResidOriR2{roinum}(isplit,ivox) = rsquared(residOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                % 
                
                %prediction of constrained model:
                % predOriResidual = voxBetas' - squeeze(voxPredFeatCoef{roinum}(isplit,ivox,:))'*voxPrfFeatSample';
                % voxPredOriR2{roinum}(isplit,ivox) = rsquared(predOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
            end
        end
        
        for isplit=1:nsplits
            imgTrials = splitImgTrials(isplit,:);
            numTrials = sum(imgTrials);
            for ivox=1:nvox(roinum)
                
                voxBetas = roiBetas{roinum}(ivox,imgTrials>0)';
                
                voxPrfSample = squeeze(prfSampleLev{roinum}(imgNum(imgTrials>0),ivox,:));
                %add constant predictor
                voxPrfSample(:,end+1) = ones;
                voxPrfFeatSample = squeeze(prfSampleLevOri{roinum}(imgNum(imgTrials>0),ivox,:,:));
                voxPrfFeatSample = reshape(voxPrfFeatSample,[],numLevels*numFeatures);
                
                %add constant predictor
                voxPrfFeatSample(:,end+1) = ones;
                
                voxFeatResidualSplit{roinum}(isplit,ivox,1:numTrials) = voxBetas' - squeeze(voxFeatCoef{roinum}(nsplits-isplit+1,ivox,:))'*voxPrfFeatSample';
                voxResidualSplit{roinum}(isplit,ivox,1:numTrials) = voxBetas' - squeeze(voxCoef{roinum}(nsplits-isplit+1,ivox,:))'*voxPrfSample';
                
                %r2 between splits
                r2split{roinum}(isplit,ivox) = rsquared(voxResidualSplit{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                r2featSplit{roinum}(isplit,ivox) = rsquared(voxFeatResidualSplit{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                aicFeatSplit{roinum}(isplit,ivox) = computeAIC(voxFeatResidualSplit{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), numTrials, size(voxFeatCoef{roinum},3));
                bicFeatSplit{roinum}(isplit,ivox) = computeBIC(voxFeatResidualSplit{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), numTrials, size(voxFeatCoef{roinum},3));
                
                %corr between splits
                pearsonRfeat{roinum}(isplit,ivox) = corr(voxBetas,(squeeze(voxFeatCoef{roinum}(nsplits-isplit+1,ivox,:))'*voxPrfFeatSample')');
                pearsonR{roinum}(isplit,ivox) = corr(voxBetas,(squeeze(voxCoef{roinum}(nsplits-isplit+1,ivox,:))'*voxPrfSample')');
                
            end
        end
        
        
    end
    
    nsd.voxResidual = voxResidual;
    nsd.voxFeatResidual = voxFeatResidual;
    nsd.voxResidualSplit = voxResidualSplit;
    nsd.voxFeatResidualSplit = voxFeatResidualSplit;
    nsd.r2 = r2;
    nsd.r2feat = r2feat;
    nsd.r2split = r2split;
    nsd.r2featSplit = r2featSplit;
    nsd.aicFeatSplit = aicFeatSplit;
    nsd.bicFeatSplit = bicFeatSplit;
    nsd.pearsonRfeat = pearsonRfeat;
    nsd.pearsonR = pearsonR;
    nsd.imgNum = imgNum;
    nsd.splitImgTrials = splitImgTrials;
    nsd.voxCoef = voxCoef;
    nsd.voxFeatCoef = voxFeatCoef;
    nsd.voxPredFeatCoef = voxPredFeatCoef;
    nsd.voxFeatPredFeatCoef = voxFeatPredFeatCoef;
    nsd.voxResidFeatCoef = voxResidFeatCoef;
    nsd.voxFeatResidFeatCoef = voxFeatResidFeatCoef;
    nsd.voxPredFeatR2 = voxPredFeatR2;
    nsd.voxFeatPredFeatR2 = voxFeatPredFeatR2;
    nsd.voxResidFeatR2 = voxResidFeatR2;
    nsd.voxFeatResidFeatR2 = voxFeatResidFeatR2;
    nsd.roiInd = roiInd;
    
    
 
    %% SAVE RESULTS
    % bandpassStr = ['_bandpass' num2str(bandMin) 'to' num2str(bandMax)];
    save(fullfile(boxfolder,['regressPrfSplit_', method, '_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']), ...
        'nsd',...
        'numLevels', 'numFeatures','rois','nvox','roiPrf','nsplits');
    toc
end

end

