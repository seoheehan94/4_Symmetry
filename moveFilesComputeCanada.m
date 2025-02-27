% cd '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli';
% nsdfolder = '/bwdata/NSDData/nsddata/experiments/nsd/';
% nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
% nsdDesign = load(nsdDesignFilename);
% allImgs = nsdDesign.sharedix; %indices of the shared 1000 images
% featurefolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/parfilter/';%to save model outputs

% mkdir pyramid/subj05
% mkdir subj06
% mkdir subj07
% % mkdir subj08
% for isub= 5:5
% 
%     allImgs = nsdDesign.subjectim(isub,nsdDesign.masterordering);%indices of all 10000 images used for this subject
%     allImgs = unique(allImgs);
%     folderName = ['pyramid/subj0', num2str(isub)];
%     for curImg = 1: length(allImgs)
%         curImgName = ['pyramid/pyrImg', num2str(allImgs(curImg)), '.mat'];
% 
%         if isfile(fullfile(curImgName))
%             fprintf('%s....\n',curImgName);
%             movefile(curImgName, folderName)
%         end
% 
%     end
    % command = 'scp -r /bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/mirfilter/ hanseohe@beluga4.computecanada.ca:/home/hanseohe/scratch/stimuli/';
    % system(command);

    % movefile([folderName, '/*'], 'pyramid/')


% end
command = ['scp -r /bwdata/NSDData/stimuli/vecLD/{images01,images02,images03,images04,images05} ' ...
           'hanseohe@beluga4.computecanada.ca:/home/hanseohe/scratch/stimuli/vecLD'];
system(command);
