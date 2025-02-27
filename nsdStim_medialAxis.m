% nsdStim.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: nsdStim()
%   by: zvi roth
%   date: 7/29/2022
%   purpose: run natural scene stimuli through the models, get filter energy responses,
%   and save output
%   creates files used by: prfSampleModel.m

% uses the steerable pyramid: https://github.com/elimerriam/stimulusVignetting

close all;
clear all;

backgroundSize = [512 512];
renderSize = [357,357];

% methods = {'contour', 'medialAxis'};
filefolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/parfilter/';%to save model outputs
fileList = dir(filefolder);
fileList = fileList(~ismember({fileList(:).name},{'.','..'}));
savefolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/medialAxisfilter/';%to save model outputs

%%
numFeatures = 1;
bandwidth = 1;
dims = backgroundSize;
numLevels = 1;

%%
iimg=0;
for imgNum=1:length(fileList)
    iimg = iimg+1
    name = fileList(imgNum).name;
    filename = ['maImg' num2str(imgNum) '.mat'];
    if ~isfile(fullfile(savefolder, filename))%if file exists already no need to remake it

        load(fullfile(filefolder, name), 'model');

        model = rmfield(model,'area');

        model.contour(model.contour ~= 0) = 1;
        model.medialAxis(model.medialAxis ~= 0) = 1;

        save(fullfile(savefolder, filename),...
            'numFeatures','bandwidth','dims','model','numLevels');
    end
end

