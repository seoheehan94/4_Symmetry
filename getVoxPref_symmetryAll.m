% getVoxPref.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: getVoxPref(1,1)
%   by: zvi roth
%   date: 7/29/2022
%   purpose: extract preferred orientation from regression weights, and sum
%   across partitions
%   uses files created by: regressPrfSplit.m
%   creates files used by: fig$$.m


function getVoxPref_symmetryAll(isub,numregions, type)

%type = {'contour', 'medialAxis', 'area'};
curPrf = ['/bwdata/NSDData/Seohee/Symmetry/prfsample_',type,'/'];

% mrQuit
close all
% clear all
global interpSz;
global backgroundSz;
global degPerPix;

numFeatures = 3;
nperms=1000;

interpSz= 714;
backgroundSz= 1024;

imgScaling = 0.5;
interpSz= 714*imgScaling;
backgroundSz= 1024*imgScaling;
degPerPix = 8.4/(714*imgScaling);
for iregion=1:numregions
    visualRegion = iregion;%V1,V2,V3,V4
    visualRegion
    load(fullfile(curPrf,['regressPrfSplit_', type,  '_v' num2str(visualRegion) '_sub' num2str(isub)  '.mat']), ...
        'nsd',...
        'numLevels', 'numFeatures','rois','nvox','roiPrf','nsplits');
    
    if length(rois)>1 %combine across ventral and dorsal ROIs
        oldNsd = nsd;
        nsd.voxResidual{1} = [];
        nsd.voxFeatResidual{1} = [];
        nsd.voxResidualSplit{1} = [];
        nsd.voxFeatResidualSplit{1} = [];
        nsd.voxFeatFstat{1} = [];
        nsd.voxFeatpvalue{1} = [];
        nsd.r2{1} = [];
        nsd.r2feat{1} = [];
        nsd.r2split{1} = [];
        nsd.r2featSplit{1} = [];
        nsd.r2featSplit_only_par{1} = [];
        nsd.r2featSplit_only_mir{1} = [];
        nsd.r2featSplit_only_tap{1} = [];
        nsd.aicFeatSplit{1} = [];
        nsd.bicFeatSplit{1} = [];
        nsd.voxCoef{1} = [];
        nsd.voxFeatCoef{1} = [];
        nsd.voxPredFeatCoef{1} = [];
        nsd.voxFeatPredFeatCoef{1} = [];
        nsd.voxResidFeatCoef{1} = [];
        nsd.voxFeatResidFeatCoef{1} = [];

        nsd.voxPredFeatR2{1} = [];
        nsd.voxFeatPredFeatR2{1} = [];
        nsd.voxResidFeatR2{1} = [];
        nsd.voxFeatResidFeatR2{1} = [];

        nsd.pearsonRfeat{1} = [];
        nsd.pearsonR{1} = [];

        nsd.roiInd{1} = [];

  
        for iroi=1:length(rois)
            nsd.voxResidual{1} = cat(2,nsd.voxResidual{1},oldNsd.voxResidual{iroi});
            nsd.voxFeatResidual{1} = cat(2,nsd.voxFeatResidual{1},oldNsd.voxFeatResidual{iroi});
            nsd.voxResidualSplit{1} = cat(2,nsd.voxResidualSplit{1},oldNsd.voxFeatResidual{iroi});
            nsd.voxFeatResidualSplit{1} = cat(2,nsd.voxFeatResidualSplit{1},oldNsd.voxFeatResidual{iroi});

            nsd.pearsonRfeat{1} = cat(2,nsd.pearsonRfeat{1},oldNsd.pearsonRfeat{iroi});
            nsd.pearsonR{1} = cat(2,nsd.pearsonR{1},oldNsd.pearsonR{iroi});
            nsd.voxFeatFstat{1} = cat(2,nsd.voxFeatFstat{1},oldNsd.voxFeatFstat{iroi});
            nsd.voxFeatpvalue{1} = cat(2,nsd.voxFeatpvalue{1},oldNsd.voxFeatpvalue{iroi});
            nsd.r2{1} = cat(2,nsd.r2{1},oldNsd.r2{iroi});
            nsd.r2feat{1} = cat(2,nsd.r2feat{1},oldNsd.r2feat{iroi});
            nsd.r2split{1} = cat(2,nsd.r2split{1},oldNsd.r2split{iroi});
            nsd.r2featSplit{1} = cat(2,nsd.r2featSplit{1},oldNsd.r2featSplit{iroi});
            nsd.r2featSplit_only_par{1} = cat(2,nsd.r2featSplit_only_par{1},oldNsd.r2featSplit_only_par{iroi});
            nsd.r2featSplit_only_mir{1} = cat(2,nsd.r2featSplit_only_mir{1},oldNsd.r2featSplit_only_mir{iroi});
            nsd.r2featSplit_only_tap{1} = cat(2,nsd.r2featSplit_only_tap{1},oldNsd.r2featSplit_only_tap{iroi});
            nsd.aicFeatSplit{1} = cat(2,nsd.aicFeatSplit{1},oldNsd.aicFeatSplit{iroi});
            nsd.bicFeatSplit{1} = cat(2,nsd.bicFeatSplit{1},oldNsd.bicFeatSplit{iroi});
            nsd.voxCoef{1} = cat(2,nsd.voxCoef{1},oldNsd.voxCoef{iroi});
            nsd.voxFeatCoef{1} = cat(2,nsd.voxFeatCoef{1},oldNsd.voxFeatCoef{iroi});
            nsd.voxPredFeatCoef{1} = cat(2,nsd.voxPredFeatCoef{1},oldNsd.voxPredFeatCoef{iroi});
            nsd.voxFeatPredFeatCoef{1} = cat(2,nsd.voxFeatPredFeatCoef{1},oldNsd.voxFeatPredFeatCoef{iroi});
            nsd.voxResidFeatCoef{1} = cat(2,nsd.voxResidFeatCoef{1},oldNsd.voxResidFeatCoef{iroi});
            nsd.voxFeatResidFeatCoef{1} = cat(2,nsd.voxFeatResidFeatCoef{1},oldNsd.voxFeatResidFeatCoef{iroi});

            nsd.voxPredFeatR2{1} = cat(2,nsd.voxPredFeatR2{1},oldNsd.voxPredFeatR2{iroi});
            nsd.voxFeatPredFeatR2{1} = cat(2,nsd.voxFeatPredFeatR2{1},oldNsd.voxFeatPredFeatR2{iroi});
            nsd.voxResidFeatR2{1} = cat(2,nsd.voxResidFeatR2{1},oldNsd.voxResidFeatR2{iroi});
            nsd.voxFeatResidFeatR2{1} = cat(2,nsd.voxFeatResidFeatR2{1},oldNsd.voxFeatResidFeatR2{iroi});

            nsd.roiInd{1} = cat(1,nsd.roiInd{1}, oldNsd.roiInd{iroi});

            %
        end
        oldPrf = roiPrf; clear roiPrf;
        roiPrf{1}.ecc=[];
        roiPrf{1}.ang=[];
        roiPrf{1}.sz=[];
        %         roiPrf{1}.exponent=[];
        %         roiPrf{1}.gain=[];
        roiPrf{1}.r2=[];
        roiPrf{1}.x=[];
        roiPrf{1}.y=[];
        for iroi=1:length(rois)
            roiPrf{1}.ecc = cat(1,roiPrf{1}.ecc,oldPrf{iroi}.ecc);
            roiPrf{1}.ang = cat(1,roiPrf{1}.ang,oldPrf{iroi}.ang);
            roiPrf{1}.sz = cat(1,roiPrf{1}.sz,oldPrf{iroi}.sz);
            %             roiPrf{1}.exponent = cat(1,roiPrf{1}.exponent,oldPrf{iroi}.exponent);
            %             roiPrf{1}.gain = cat(1,roiPrf{1}.gain,oldPrf{iroi}.gain);
            roiPrf{1}.r2 = cat(1,roiPrf{1}.r2,oldPrf{iroi}.r2);
            roiPrf{1}.x = cat(1,roiPrf{1}.x,oldPrf{iroi}.x);
            roiPrf{1}.y = cat(1,roiPrf{1}.y,oldPrf{iroi}.y);
        end
        rois = 1;
    end

    %% AVERAGE SPLITS
    nsd.voxCoef{1}(nsplits+1,:,:) = mean(nsd.voxCoef{1},1);
    nsd.voxFeatCoef{1}(nsplits+1,:,:) = mean(nsd.voxFeatCoef{1},1);
    nsd.voxPredFeatCoef{1}(nsplits+1,:,:) = mean(nsd.voxPredFeatCoef{1},1);
    nsd.voxFeatPredFeatCoef{1}(nsplits+1,:,:) = mean(nsd.voxFeatPredFeatCoef{1},1);
    nsd.voxResidFeatCoef{1}(nsplits+1,:,:) = mean(nsd.voxResidFeatCoef{1},1);
    nsd.voxFeatResidFeatCoef{1}(nsplits+1,:,:) = mean(nsd.voxFeatResidFeatCoef{1},1);

    nsd.pearsonRfeat{1}(nsplits+1,:) = mean(nsd.pearsonRfeat{1},1);
    nsd.pearsonR{1}(nsplits+1,:) = mean(nsd.pearsonR{1},1);
    nsd.voxFeatFstat{1} = mean(nsd.voxFeatFstat{1},3);
    nsd.voxFeatpvalue{1}= mean(nsd.voxFeatpvalue{1},3);
    nsd.voxFeatFstat{1}(nsplits+1,:) = mean(nsd.voxFeatFstat{1},1);
    nsd.voxFeatpvalue{1}(nsplits+1,:)= mean(nsd.voxFeatpvalue{1},1);
    nsd.r2{1}(nsplits+1,:) = mean(nsd.r2{1},1);
    nsd.r2feat{1}(nsplits+1,:) = mean(nsd.r2feat{1},1);
    nsd.r2split{1}(nsplits+1,:) = mean(nsd.r2split{1},1);
    nsd.r2featSplit{1}(nsplits+1,:) = mean(nsd.r2featSplit{1},1);
    nsd.r2featSplit_only_par{1}(nsplits+1,:) = mean(nsd.r2featSplit_only_par{1},1);
    nsd.r2featSplit_only_mir{1}(nsplits+1,:) = mean(nsd.r2featSplit_only_mir{1},1);
    nsd.r2featSplit_only_tap{1}(nsplits+1,:) = mean(nsd.r2featSplit_only_tap{1},1);
    nsd.aicFeatSplit{1}(nsplits+1,:) = mean(nsd.aicFeatSplit{1},1);
    nsd.bicFeatSplit{1}(nsplits+1,:) = mean(nsd.bicFeatSplit{1},1);


    nsplits = nsplits+1;
  

    %% COMPUTE PREFERRED ORIENTATION - CIRCULAR CENTER OF MASS
    % clear fullPrefOri residPrefOri residOriPrefOri predCoefOri predOriCoefOri oriDeviation vertDeviation cardDeviation
    % clear fullOriModul fullPrefAmp fullAntiAmp
    % for  iroi=1:length(rois)%rois=1
    %     for isplit=1:nsplits
    %         %full orientation model
    %         fullCoef = squeeze(nsd.voxOriCoef{iroi}(isplit,:,1:end-1));
    %         fullPrefOri{iroi}(isplit,:) = gratingPrefOri(fullCoef,gratings.modelOriEnergy);
    % 
    %         %full model on residuals of constrained model
    %         fullCoef = squeeze(nsd.voxResidOriCoef{iroi}(isplit,:,1:end-1));
    %         residPrefOri{iroi}(isplit,:) = gratingPrefOri(fullCoef,gratings.modelOriEnergy);
    % 
    %         %full model on residual of full model
    %         fullCoef = squeeze(nsd.voxOriResidOriCoef{iroi}(isplit,:,1:end-1));
    %         residOriPrefOri{iroi}(isplit,:) = gratingPrefOri(fullCoef,gratings.modelOriEnergy);
    % 
    %         %full model on constrained prediction
    %         fullCoef = squeeze(nsd.voxPredOriCoef{iroi}(isplit,:,1:end-1));
    %         predCoefOri{iroi}(isplit,:) = gratingPrefOri(fullCoef,gratings.modelOriEnergy);
    % 
    %         %full model on full model prediction
    %         fullCoef = squeeze(nsd.voxOriPredOriCoef{iroi}(isplit,:,1:end-1));
    %         predOriCoefOri{iroi}(isplit,:) = gratingPrefOri(fullCoef,gratings.modelOriEnergy);
    % 
    %     end
    %     %cross-validated measure of orientation modulation
    %     for isplit=1:2
    %         fullCoef = squeeze(nsd.voxOriCoef{iroi}(isplit,:,1:end-1));
    %         [fullOriModul{iroi}(isplit,:), fullPrefAmp{iroi}(isplit,:), fullAntiAmp{iroi}(isplit,:)] = gratingOriModulation(fullPrefOri{iroi}(3-isplit,:),fullCoef,gratings.modelOriEnergy);
    %     end
    %     fullOriModul{iroi}(nsplits,:) = mean(fullOriModul{iroi}(1:2,:),1);
    %     fullPrefAmp{iroi}(nsplits,:) = mean(fullPrefAmp{iroi}(1:2,:),1);
    %     fullAntiAmp{iroi}(nsplits,:) = mean(fullAntiAmp{iroi}(1:2,:),1);
    % 
    %     %%%%%%%%%%%%% SYNTH
    %     %full orientation model
    %     %         synthFullCoef = squeeze(synth.voxOriCoef{iroi}(:,1:end-1));
    %     %         synthFullPrefOri{iroi} = gratingPrefOri(synthFullCoef,gratings.modelOriEnergy);
    % end


    %% get symmetry preference 
    % clear fullPrefFeat fullCoef
    % for  iroi=1:length(rois)%rois=1
    %     for isplit=1:nsplits
    %         fullFstat = squeeze(nsd.voxFeatFstat{iroi}(isplit,:));
    %         numvox = size(fullFstat,2);
    % 
    %         % coefMat = reshape(fullCoef,numvox,numLevels,numFeatures);%vox x levels x orientations
    % 
    %         fullPrefFeat{iroi}(isplit,:) = fullFstat;
    %     end
    % end
    % 

   
    %%
    %save preferred orientation and level for this ROI
    iroi=1;
    allRoiPrf{iregion} = roiPrf{iroi};%iroi=1
    roiFstat{iregion} = nsd.voxFeatFstat{iroi};
    roiFstat_p{iregion} = nsd.voxFeatpvalue{iroi};

    roiInd{iregion} = nsd.roiInd{iroi};
    roiNsdFeatR2{iregion} = nsd.r2featSplit{iroi};
    roiNsdFeatR2_par{iregion} = nsd.r2featSplit_only_par{iroi};
    roiNsdFeatR2_mir{iregion} = nsd.r2featSplit_only_mir{iroi};
    roiNsdFeatR2_tap{iregion} = nsd.r2featSplit_only_tap{iroi};
    allaicFeatSplit{iregion} = nsd.aicFeatSplit{iroi};
    allbicFeatSplit{iregion} = nsd.bicFeatSplit{iroi};

end


%save all ROIs to create overlay
roifolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
visualRoisFile = fullfile(roifolder,'roi/prf-visualrois.nii.gz');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4'};
combinedRoiNames = {'V1','V2','V3','hV4'};

save([curPrf, 'voxModelPref_', type, '_sub', num2str(isub), '.mat'],'allRoiPrf',...
    'roiFstat', 'roiFstat_p', 'roiNsdFeatR2', 'roiNsdFeatR2_par', 'roiNsdFeatR2_mir', 'roiNsdFeatR2_tap',...
    'allaicFeatSplit', 'allbicFeatSplit',...
    'visRoiData','roiNames','combinedRoiNames','roiInd','nsplits');


