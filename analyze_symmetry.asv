%analyze_orientation.m
% get voxel preference from model weights and create brainVolume
%   uses files created by: regressPrfSplit.m
%   creates files used by:
clear all;
type = 'par'; %'par', 'mir', 'taper'
savefolder = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/brainVolume_', type];
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4','OPA','PPA','RSC'};
combinedRoiNames = {'V1','V2','V3','hV4','OPA','PPA','RSC'};
methods = {'contour', 'medialAxis', 'area'};
curPrf = ['/bwdata/NSDData/Seohee/Symmetry/prfsample_' type, '/'];

%% 1.  mean R2, AIC, BIC %%
totalR2FeatSplit = struct;
totalaicFeatSplit = struct;
totalbicFeatSplit = struct;
for method = 1:3
    totalR2FeatSplit.(methods{method}) = {};
    totalaicFeatSplit.(methods{method}) = {};
    totalbicFeatSplit.(methods{method}) = {};
    for isub = 1:8
        fprintf('isub:%d. method:%d. ...\n',isub,method);
        load([curPrf 'voxModelPref_' methods{method} '_sub' num2str(isub) '.mat']);

        %% total values of R2, aic, bic
        totalR2FeatSplit.(methods{method}){end+1} = roiNsdFeatR2;
        totalaicFeatSplit.(methods{method}){end+1} = allaicFeatSplit;
        totalbicFeatSplit.(methods{method}){end+1} = allbicFeatSplit;

    end
end

%% R2
fieldsCon = fieldnames(totalR2FeatSplit);
allroiR2FeatSplit=[];
V4R2FeatSplit = [];
OPAR2FeatSplit=[];
PPAR2FeatSplit=[];
RSCR2FeatSplit=[];
for i = 1:numel(fieldsCon)
    curRoiR2FeatSplit = [];
    curV4R2FeatSplit = [];
    curOPAR2FeatSplit = [];
    curPPAR2FeatSplit = [];
    curRSCR2FeatSplit = [];

    for sub = 1:numel(totalR2FeatSplit.(fieldsCon{i}))
        for roi = 1:size(roiNsdFeatR2,2)
            curRoiR2FeatSplit = [curRoiR2FeatSplit, totalR2FeatSplit.(fieldsCon{i}){sub}{roi}(3,:)];
        end
        curV4R2FeatSplit = [curV4R2FeatSplit, totalR2FeatSplit.(fieldsCon{i}){sub}{4}(3,:)];
        curOPAR2FeatSplit = [curOPAR2FeatSplit, totalR2FeatSplit.(fieldsCon{i}){sub}{5}(3,:)];
        curPPAR2FeatSplit = [curPPAR2FeatSplit, totalR2FeatSplit.(fieldsCon{i}){sub}{6}(3,:)];
        curRSCR2FeatSplit = [curRSCR2FeatSplit, totalR2FeatSplit.(fieldsCon{i}){sub}{7}(3,:)];
    end
    writematrix(curRoiR2FeatSplit', ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/analyses_', type, '/allroiR2', '_', (fieldsCon{i}), '.csv']);
    writematrix(curV4R2FeatSplit', ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/analyses_', type, '/V4R2_', (fieldsCon{i}), '.csv']);
    writematrix(curOPAR2FeatSplit', ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/analyses_', type, '/OPAR2_', (fieldsCon{i}), '.csv']);
    writematrix(curPPAR2FeatSplit', ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/analyses_', type, '/PPAR2_', (fieldsCon{i}), '.csv']);
    writematrix(curRSCR2FeatSplit', ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/analyses_', type, '/RSCR2_', (fieldsCon{i}), '.csv']);

    allroiR2FeatSplit(i) = mean(curRoiR2FeatSplit, 'omitnan');
    V4R2FeatSplit(i) = mean(curV4R2FeatSplit,'omitnan');
    OPAR2FeatSplit(i) = mean(curOPAR2FeatSplit,'omitnan');
    PPAR2FeatSplit(i) = mean(curPPAR2FeatSplit,'omitnan');
    RSCR2FeatSplit(i) = mean(curRSCR2FeatSplit,'omitnan');
end

save(fullfile(['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/analyses_', type, '/', 'meanR2.mat']), "allroiR2FeatSplit", "V4R2FeatSplit", "OPAR2FeatSplit", "PPAR2FeatSplit", "RSCR2FeatSplit");


%% AIC/BIC
fieldsCon = fieldnames(totalaicFeatSplit);
for i = 1:numel(fieldsCon)
    curRoiaicFeatSplit = [];
    curRoibicFeatSplit = [];
    curV4aicFeatSplit = [];
    curV4bicFeatSplit = [];
    curOPAaicFeatSplit = [];
    curOPAbicFeatSplit = [];
    curPPAaicFeatSplit = [];
    curPPAbicFeatSplit = [];
    curRSCaicFeatSplit = [];
    curRSCbicFeatSplit = [];
    for sub = 1:numel(totalR2FeatSplit.(fieldsCon{i}))
        for roi = 1:size(roiNsdFeatR2,2)
            curRoiaicFeatSplit = [curRoiaicFeatSplit, totalaicFeatSplit.(fieldsCon{i}){sub}{roi}(3,:)];
            curRoibicFeatSplit = [curRoibicFeatSplit, totalbicFeatSplit.(fieldsCon{i}){sub}{roi}(3,:)];
        end

        curV4aicFeatSplit = [curV4aicFeatSplit, totalaicFeatSplit.(fieldsCon{i}){sub}{4}(3,:)];
        curV4bicFeatSplit = [curV4bicFeatSplit, totalbicFeatSplit.(fieldsCon{i}){sub}{4}(3,:)];
        curOPAaicFeatSplit = [curOPAaicFeatSplit, totalaicFeatSplit.(fieldsCon{i}){sub}{5}(3,:)];
        curOPAbicFeatSplit = [curOPAbicFeatSplit, totalbicFeatSplit.(fieldsCon{i}){sub}{5}(3,:)];
        curPPAaicFeatSplit = [curPPAaicFeatSplit, totalaicFeatSplit.(fieldsCon{i}){sub}{6}(3,:)];
        curPPAbicFeatSplit = [curPPAbicFeatSplit, totalbicFeatSplit.(fieldsCon{i}){sub}{6}(3,:)];
        curRSCaicFeatSplit = [curRSCaicFeatSplit, totalaicFeatSplit.(fieldsCon{i}){sub}{7}(3,:)];
        curRSCbicFeatSplit = [curRSCbicFeatSplit, totalbicFeatSplit.(fieldsCon{i}){sub}{7}(3,:)];
    end
    % writematrix(curRoiaicFeatSplit, ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/allroiR2', imgType{curimgtype}, '_', (fieldsCon{i}), '.csv']);
    % writematrix(curRoiaicFeatSplit, ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/MaxMin/V1R2', imgType{curimgtype}, '_', (fieldsCon{i}), '.csv']);

    allroiaicFeatSplit(i) = mean(curRoiaicFeatSplit,"omitnan");
    allroibicFeatSplit(i) = mean(curRoibicFeatSplit,"omitnan");
    V4aicFeatSplit(i) = mean(curV4aicFeatSplit,"omitnan");
    V4bicFeatSplit(i) = mean(curV4bicFeatSplit,"omitnan");
    OPAaicFeatSplit(i) = mean(curOPAaicFeatSplit,"omitnan");
    OPAbicFeatSplit(i) = mean(curOPAbicFeatSplit,"omitnan");
    PPAaicFeatSplit(i) = mean(curPPAaicFeatSplit,"omitnan");
    PPAbicFeatSplit(i) = mean(curPPAbicFeatSplit,"omitnan");
    RSCaicFeatSplit(i) = mean(curRSCaicFeatSplit,"omitnan");
    RSCbicFeatSplit(i) = mean(curRSCbicFeatSplit,"omitnan");
end

save(fullfile(['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/analyses_', type, '/', 'meanAICBIC.mat']), "allroiaicFeatSplit", "allroibicFeatSplit", ...
    "V4aicFeatSplit","V4bicFeatSplit", ...
    "OPAaicFeatSplit", "OPAaicFeatSplit",...
    "PPAaicFeatSplit", "PPAaicFeatSplit", ...
    "RSCaicFeatSplit", "RSCbicFeatSplit");


%% Count positive R2

positiveCount_sub = struct;
totalCount_sub = struct;
for method = 1:length(methods)
    positiveCount_sub.(methods{method}) = {};
    NACount_sub.(methods{method}) = {};
    totalCount_sub.(methods{method}) = {};
    for isub = 1:8
         for roi = 1:size(roiNsdFeatR2,2)
            curR2 =  totalR2FeatSplit.(methods{method}){isub}{roi}(3,:);
            positiveCount_sub.(methods{method}){isub}(roi) = sum(curR2(:)>=0);
            NACount_sub.(methods{method}){isub}(roi) = sum(isnan(curR2(:)));
            totalCount_sub.(methods{method}){isub}(roi) = size(curR2,2);
         end
    end
end

positiveCount = struct;
for method = 1:length(methods)
    positiveCount.(methods{method}) = [];
    for roi = 1:7
        curRoiPositiveSum =0;
        curRoiTotalSum =0;
        curRoiNASum =0;
        for isub = 1:8
            curRoiPositiveSum = curRoiPositiveSum+positiveCount_sub.(methods{method}){isub}(roi);
            curRoiNASum = curRoiNASum+NACount_sub.(methods{method}){isub}(roi);
            curRoiTotalSum = curRoiTotalSum+totalCount_sub.(methods{method}){isub}(roi);
        end
        positiveCount.(methods{method})(1,roi) = curRoiPositiveSum;
        positiveCount.(methods{method})(2,roi) = curRoiTotalSum;
        positiveCount.(methods{method})(3,roi) = curRoiNASum;
    end
end

% save(fullfile(['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/analyses_', type, '/', 'positiveR2Count.mat']), "positiveCount")
%% make a brain volume

for method = 1:3

    for isub = 1:8
        fprintf('isub:%d. con:%d. ...\n',isub,method);
        load([curPrf 'voxModelPref_' methods{method} '_sub' num2str(isub) '.mat']);
        % save all ROIs to create overlay
        roifolder = ['/bwdata/NSDData/nsddata/ppdata/subj0' num2str(isub) '/func1pt8mm/'];
        visualRoisFile = fullfile(roifolder,'roi/prf-visualrois.nii.gz');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
        visRoiData = niftiread(visualRoisFile);

        placesRoisFile = fullfile(roifolder,'roi/floc-places.nii.gz'); %OPA, PPA, RSC
        placeRoiData = niftiread(placesRoisFile);
        visRoiData(placeRoiData == 1) = 8;
        visRoiData(placeRoiData == 2) = 9;
        visRoiData(placeRoiData == 3) = 10;

        ourBrain = visRoiData;
        ourBrain(ourBrain == 2) = 1;
        ourBrain(ourBrain == 3 | ourBrain == 4) = 2;
        ourBrain(ourBrain == 5 | ourBrain == 6) = 3;
        ourBrain(ourBrain == 7) = 4;
        ourBrain(ourBrain == 8) = 5;
        ourBrain(ourBrain == 9) = 6;
        ourBrain(ourBrain == 10) = 7;


        % make a brain volume
        newBrain = ourBrain;
        newBrain(newBrain <= 0) = NaN;
        newBrain(newBrain > 0) = 0;
        for visualRegion = 1:7
            curOurBrain = ourBrain;
            % if visualRegion == 2
            %     curOurBrain(visRoiData == 3 | visRoiData == 4) = 2;
            % elseif visualRegion == 3
            %     curOurBrain(visRoiData == 5 | visRoiData == 6) = 3;
            % end
            newBrain(curOurBrain == visualRegion) = roiFeat{visualRegion}(3,:);
        end
        %
        % for visualRegion = 1:4
        %     curOurBrain = ourBrain;
        %     % if visualRegion == 2
        %     %     curOurBrain(visRoiData == 3 | visRoiData == 4) = 2;
        %     % elseif visualRegion == 3
        %     %     curOurBrain(visRoiData == 5 | visRoiData == 6) = 3;
        %     % end
        %     curNewBrain = curOurBrain;
        %     curNewBrain(curOurBrain ~= visualRegion) = -1;
        %     thisfield = combinedRoiNames{visualRegion};
        %
        %     curNewBrain(curOurBrain == visualRegion) = roiOri.(thisfield);
        %
        %     newBrainbyROI(:,:,:,visualRegion) =curNewBrain;
        % end

        % save(fullfile(savefolder, [condition{con}, 'Brain_sub', num2str(isub), '.mat']), 'newBrain');
        % save(fullfile(savefolder, [condition{con}, 'BrainbyROI_sub', num2str(isub), '.mat']), 'newBrainbyROI');

        % R2
        r2Brain = ourBrain;
        r2Brain(ourBrain <= 0) = NaN;
        r2Brain(ourBrain > 0) = 0;

        for visualRegion = 1:7
            r2Brain(ourBrain == visualRegion) = roiNsdFeatR2{visualRegion}(3,:);
        end

        %
        % % %% save afni file
        % % % size(newBrain)
        % % %3dinfo betas_session01.nii.gz
        % % %3dcalc -a betas_session01.nii.gz[1] -expr a -prefix oneBeta
        % % % command = ['3dinfo betas_session01_sub', num2str(isub), '.nii.gz'];
        % % % system(command);
        % % % command = ['3dcalc -a betas_session01_sub', num2str(isub), '.nii.gz[1] -expr a -prefix oneBeta_sub', num2str(isub)];
        % % % system(command);
        % %
        % % currOneBeta = ['oneBeta_sub', num2str(isub), '+orig'];
        % % [err,V,Info] = BrikLoad(currOneBeta);
        % %
        % % Info.RootName = ['curvBrain_sub', num2str(isub), '+orig'];
        % % opt.Prefix = ['curvBrain_sub', num2str(isub)];
        % % WriteBrik(newBrain,Info,opt);
        % % Info.RootName = ['curvBrainbyROI_sub', num2str(isub), '+orig'];
        % % opt.Prefix = ['curvBrainbyROI_sub', num2str(isub)];
        % % WriteBrik(newBrainbyROI,Info,opt);
        %
        %% save nifti
        % cd('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Curvature/brainVolume_regress');
        % load(['oriBrain_sub', num2str(isub), '.mat']);
        info_old = niftiinfo(['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/brainVolume/betas_session01_sub', num2str(isub),'.nii.gz']);

        niftiwrite(newBrain,[savefolder, '/', methods{method}, 'Brain_sub', num2str(isub),'.nii']);
        info_new = niftiinfo([savefolder, '/', methods{method}, 'Brain_sub', num2str(isub),'.nii']);
        info_new.PixelDimensions = [1.8, 1.8, 1.8];
        info_new.TransformName = info_old.TransformName;
        info_new.SpatialDimension = info_old.SpatialDimension;
        info_new.Transform = info_old.Transform;
        info_new.Qfactor = info_old.Qfactor;
        info_new.AuxiliaryFile = info_old.AuxiliaryFile;
        info_new.raw.pixdim = info_old.raw.pixdim;
        info_new.raw.aux_file = info_old.raw.aux_file;
        info_new.raw.sform_code = info_old.raw.sform_code;
        info_new.raw.srow_x = info_old.raw.srow_x;
        info_new.raw.srow_y = info_old.raw.srow_y;
        info_new.raw.srow_z = info_old.raw.srow_z;
        niftiwrite(newBrain,[savefolder, '/', methods{method}, 'Brain_sub', num2str(isub),'.nii'], info_new);

        % niftiwrite(newBrainbyROI,[savefolder, '/', condition{con}, 'BrainbyROI_sub', num2str(isub),'.nii']);
        % info_new = niftiinfo([savefolder, '/', condition{con}, 'BrainbyROI_sub', num2str(isub),'.nii']);
        % info_new.PixelDimensions = info_old.PixelDimensions;
        % info_new.TransformName = info_old.TransformName;
        % info_new.SpatialDimension = info_old.SpatialDimension;
        % info_new.Transform = info_old.Transform;
        % info_new.Qfactor = info_old.Qfactor;
        % info_new.AuxiliaryFile = info_old.AuxiliaryFile;
        % info_new.raw.pixdim = info_old.raw.pixdim;
        % info_new.raw.aux_file = info_old.raw.aux_file;
        % info_new.raw.sform_code = info_old.raw.sform_code;
        % info_new.raw.srow_x = info_old.raw.srow_x;
        % info_new.raw.srow_y = info_old.raw.srow_y;
        % info_new.raw.srow_z = info_old.raw.srow_z;
        % niftiwrite(newBrainbyROI,[savefolder, '/', condition{con}, 'BrainbyROI_sub', num2str(isub),'.nii'],info_new);


        niftiwrite(r2Brain,[savefolder, '/', methods{method}, 'BrainR2_sub', num2str(isub),'.nii']);
        info_new = niftiinfo([savefolder, '/', methods{method}, 'BrainR2_sub', num2str(isub),'.nii']);
        info_new.PixelDimensions = [1.8, 1.8, 1.8];
        info_new.TransformName = info_old.TransformName;
        info_new.SpatialDimension = info_old.SpatialDimension;
        info_new.Transform = info_old.Transform;
        info_new.Qfactor = info_old.Qfactor;
        info_new.AuxiliaryFile = info_old.AuxiliaryFile;
        info_new.raw.pixdim = info_old.raw.pixdim;
        info_new.raw.aux_file = info_old.raw.aux_file;
        info_new.raw.sform_code = info_old.raw.sform_code;
        info_new.raw.srow_x = info_old.raw.srow_x;
        info_new.raw.srow_y = info_old.raw.srow_y;
        info_new.raw.srow_z = info_old.raw.srow_z;
        niftiwrite(r2Brain,[savefolder, '/', methods{method}, 'BrainR2_sub', num2str(isub),'.nii'],info_new);



    end
end


% Load CSV files
% old_R2 = readtable('allroiR2old.csv');
% control_R2 = readtable('allroiR2control.csv');
% new_R2 = readtable('allroiR2ori.csv');
%
% % Convert table to arrays if needed
% old_R2 = old_R2{1, :}';  % Convert to column vector
% control_R2 = control_R2{1, :}';  % Convert to column vector
% new_R2 = new_R2{1, :}';  % Convert to column vector
% data = [old_R2, control_R2, new_R2];
%
%
% writematrix(old_R2, '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/allroiR2old.csv');
% writematrix(control_R2, '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/allroiR2control.csv');
% writematrix(new_R2, '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/allroiR2ori.csv');
%%
