% getGroupMean.m
% For the group map, we computed the circular mean across subjects for each vertex in V1–V4, 
% weighted by the full model R2 values.
clear all;
type = 'contour'; %'contour', 'medialAxis', 'area'

filedir = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/surfaceData_' type, '/'];
methods = {'contour', 'medialAxis', 'area'};

addpath(genpath('/usr/local/freesurfer/7.4.1/matlab'));
[~,M,mr] = load_mgh(['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/surfaceData_' type, '/sub1/' type, 'Brain_sub1_lh_fsaverage.mgh']);

% cutoff
cutoff_lh = 0;
cutoff_rh = 0;

% number of voxels survived
cutNumVoxel_lh = [];
cutNumVoxel_rh = [];

for curcond = 1:3
    %% get all subjects volume 
    % what to do with -1 and 0?
    for isub=1:8
        vol_lh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', type, 'Brain_sub', num2str(isub), '_lh_fsaverage.mgh']);
        vol_rh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', type, 'Brain_sub', num2str(isub), '_rh_fsaverage.mgh']);
    end
    

    %% get all subjects full model R2 values
    % what to do with -1 and 0?
    for isub=1:8
        R2_lh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', type, 'BrainR2_sub', num2str(isub), '_lh_fsaverage.mgh']);
        R2_rh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', type, 'BrainR2_sub', num2str(isub), '_rh_fsaverage.mgh']);
    end

 %% get all subjects full model p values
    % what to do with -1 and 0?
    for isub=1:8
        p_lh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', type, 'Brain_p_sub', num2str(isub), '_lh_fsaverage.mgh']);
        p_rh(:,isub) = load_mgh([filedir, 'sub', num2str(isub), '/', type, 'Brain_p_sub', num2str(isub), '_rh_fsaverage.mgh']);
    end
    %% weighted mean
    symPref_lh = [];
    symPref_rh = [];
    symPref_p_lh = [];
    symPref_p_rh = [];
    
    % Set Negative R² Values to Zero
    R2_lh_nonnegative = R2_lh;
    R2_rh_nonnegative = R2_rh;
    R2_lh_nonnegative(R2_lh_nonnegative < 0) = 0;
    R2_rh_nonnegative(R2_rh_nonnegative < 0) = 0;

    for ivox=1:size(vol_lh,1)
        weights = R2_lh_nonnegative(ivox, :).';
        curVol_lh = vol_lh(ivox, :).';
        curp_lh = p_lh(ivox, :).';
        % Create a logical mask where both curVol_lh and weights are non-NaN
        validMask = ~isnan(curVol_lh) & ~isnan(weights);
        validMask_p = ~isnan(curp_lh) & ~isnan(weights);
        % Apply the mask to filter out NaN values
        filteredVol = curVol_lh(validMask);
        filteredWeights = weights(validMask);

        filteredp = curp_lh(validMask_p);
        filteredWeights_p = weights(validMask_p);

        if all(isnan(weights))
            symPref_lh(ivox) = NaN;
            symPref_p_lh(ivox) = NaN;
        else
            symPref_lh(ivox) = mean(filteredVol, 'Weights', filteredWeights);
            symPref_p_lh(ivox) = mean(filteredp, 'Weights', filteredWeights_p);
        end
    end
    for ivox=1:size(vol_rh,1)
        weights = R2_rh_nonnegative(ivox, :).';
        curVol_rh = vol_rh(ivox, :).';
        curp_rh = p_rh(ivox, :).';
        % Create a logical mask where both curVol_lh and weights are non-NaN
        validMask = ~isnan(curVol_rh) & ~isnan(weights);
        validMask_p = ~isnan(curp_rh) & ~isnan(weights);
        % Apply the mask to filter out NaN values
        filteredVol = curVol_rh(validMask);
        filteredWeights = weights(validMask);

        filteredp = curp_rh(validMask_p);
        filteredWeights_p = weights(validMask_p);


        if all(isnan(weights))
            symPref_rh(ivox) = NaN;
            symPref_p_rh(ivox) = NaN;
        else
            symPref_rh(ivox) = mean(filteredVol, 'Weights', filteredWeights);
            symPref_p_rh(ivox) = mean(filteredp, 'Weights', filteredWeights_p);
        end
    end
   
   

    %% only R2 > 0 
    meanR2_lh = mean(R2_lh,2, "omitnan");
    meanR2_rh = mean(R2_rh,2, "omitnan");
    % cutoff_lh = prctile(meanR2_lh, 50);
    % cutoff_rh = prctile(meanR2_rh, 50);
    R2cutoff_symPref_lh = symPref_lh;
    R2cutoff_symPref_lh(meanR2_lh < cutoff_lh) = NaN;
    R2cutoff_symPref_rh = symPref_rh;
    R2cutoff_symPref_rh(meanR2_rh < cutoff_rh) = NaN;

    R2cutoff_symPref_p_lh = symPref_p_lh;
    R2cutoff_symPref_p_lh(meanR2_lh < cutoff_lh) = NaN;
    R2cutoff_symPref_p_rh = symPref_p_rh;
    R2cutoff_symPref_p_rh(meanR2_rh < cutoff_rh) = NaN;
    % R2cutoff_symPref_lh(R2cutoff_symPref_lh>100)= NaN;
    % R2cutoff_symPref_rh(R2cutoff_symPref_rh>100)= NaN;
    % R2cutoff_symPref_lh(R2cutoff_symPref_lh<-100)= NaN;
    % R2cutoff_symPref_rh(R2cutoff_symPref_rh<-100)= NaN;
    % meanR2_lh_save=(meanR2_lh(~isnan(meanR2_lh)));
    % meanR2_rh_save=(meanR2_rh(~isnan(meanR2_rh)));

    % writematrix(meanR2_lh_save, ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/subMeanR2_lh', condition{curcond}, '.csv']);
    % writematrix(meanR2_rh_save, ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/analyses/subMeanR_rh', condition{curcond}, '.csv']);

    %% save
    fileName_lh = [filedir, type, 'Brain_groupmean_lh_fsaverage.mgh'];
    fileName_rh = [filedir, type, 'Brain_groupmean_rh_fsaverage.mgh'];
    save_mgh(symPref_lh, fileName_lh, M,mr);
    save_mgh(symPref_rh, fileName_rh, M,mr);

    fileName_lh = [filedir, type, 'Brain_groupmean_R2cut_lh_fsaverage.mgh'];
    fileName_rh = [filedir, type, 'Brain_groupmean_R2cut_rh_fsaverage.mgh'];
    save_mgh(R2cutoff_symPref_lh, fileName_lh, M,mr);
    save_mgh(R2cutoff_symPref_rh, fileName_rh, M,mr);

    fileName_lh = [filedir, type, 'Brain_p_groupmean_lh_fsaverage.mgh'];
    fileName_rh = [filedir, type, 'Brain_p_groupmean_rh_fsaverage.mgh'];
    save_mgh(symPref_p_lh, fileName_lh, M,mr);
    save_mgh(symPref_p_rh, fileName_rh, M,mr);

    fileName_lh = [filedir, type, 'Brain_p_groupmean_R2cut_lh_fsaverage.mgh'];
    fileName_rh = [filedir, type, 'Brain_p_groupmean_R2cut_rh_fsaverage.mgh'];
    save_mgh(R2cutoff_symPref_p_lh, fileName_lh, M,mr);
    save_mgh(R2cutoff_symPref_p_rh, fileName_rh, M,mr);
end




% plot number of voxels

% figure;
% plot(numVoxel_lh,'-o', 'Color','#0072BD','LineWidth',2);
% hold on;
% plot(numVoxel_rh,'-o', 'Color','#F35872','LineWidth',2);
% xticks(1:3);
% xticklabels({'photograph-filter','line drawing-filter','contour'});
% xlabel('Methods');
% legend('left','right','FontSize',15,'FontName','Helvetica');
% ax = gca;
% ax.YAxis.FontSize = 15;
% ax.YAxis.FontName = 'Helvetica';
% title('Number of voxels above R2 0.02', 'FontSize',20,'FontName','Helvetica');
% box off;
% legend boxoff



% %% get R2 difference
% clear all;
% 
% filedir = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/';
% 
% addpath(genpath('/usr/local/freesurfer/7.4.1/matlab'));
% [~,M,mr] = load_mgh('/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/surfaceData/sub1/oldBrain_sub1_lh_fsaverage.mgh');
% 
% for isub=1:8
%     for curcond = 1:3
%         R2_lh(:,curcond) = load_mgh([filedir, 'sub', num2str(isub), '/', methods{curcond}, 'BrainR2_sub', num2str(isub), '_lh_fsaverage.mgh']);
%         R2_rh(:,curcond) = load_mgh([filedir, 'sub', num2str(isub), '/', methods{curcond}, 'BrainR2_sub', num2str(isub), '_rh_fsaverage.mgh']); 
%     end
% 
%     %old_ori
%     R2diff_lh(:,isub) = R2_lh(:,3)- R2_lh(:,1);
%     R2diff_rh(:,isub) = R2_rh(:,3)- R2_rh(:,1);
% 
% end
% 
% meanR2diff_lh = mean(R2diff_lh,2, "omitnan");
% meanR2diff_rh = mean(R2diff_rh,2, "omitnan");
% 
% fileName_lh = [filedir, 'R2diff_oriold_lh_fsaverage.mgh'];
% fileName_rh = [filedir, 'R2diff_oriold_rh_fsaverage.mgh'];
% 
% % save_mgh(meanR2diff_lh, fileName_lh, M,mr);
% % save_mgh(meanR2diff_rh, fileName_rh, M,mr);
% 
% 
% %% save by ROI
% % Path to your .label file
% visualregions = {'V1', 'V2', 'V3', 'V4'};
% hemis = {'lh', 'rh'};
% numVertices = 163842;
% for visualregion = 1:4
%     for hemi = 1:2
%         labelFilePath = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Orientation/visualRegion/labels/';
%         curfile = [labelFilePath, hemis{hemi}, '.', visualregions{visualregion}, '.label'];
% 
%         % Read the .label file
%         labelData = read_label(curfile);
% 
%         % Initialize the overlay with NaN (no value) for all vertices
%         overlay = nan(numVertices, 1);
% 
%         % Assign a value (e.g., 1.0) to vertices in the label file
%         overlay(labelData.vertexIndex) = labelData.values;
% 
%         if hemi == 1
%             curmeanR2diff = meanR2diff_lh;
%         else 
%             curmeanR2diff = meanR2diff_rh;
%         end
% 
%         curmeanR2diff(isnan(overlay))=NaN;
% 
%         fileName = [filedir, 'R2diff_oriold_', visualregions{visualregion}, '_', hemis{hemi}, '_fsaverage.mgh'];
%         save_mgh(curmeanR2diff, fileName, M,mr);
% 
%     end
% end