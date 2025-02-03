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

prompt = "Type of symmetry?";
symmetryType = input(prompt,"s");

backgroundSize = [512 512];
renderSize = [357,357];

methods = {'contour', 'medialAxis', 'area'};
% cd '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli';
%pyramidfolder = '/misc/data18/rothzn/nsd/stimuli/pyramid/';%to save model outputs
switch (lower(symmetryType))
    case 'parallelism'
        whichtype = 'par';
    case 'separation'
        whichtype = 'sep';
    case 'mirror'
        whichtype = 'mir';
    case 'taper'
        whichtype = 'tap';
end
savefolder = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/',whichtype,'filter/'];%to save model outputs

%%

numFeatures = 1;
bandwidth = 1;
dims = backgroundSize;
numLevels = 1;

%%
stimfolder = '/bwdata/NSDData/stimuli/';
stimfilename = fullfile(stimfolder,'nsd_stimuli.hdf5');%[3 1360 714 220]
stiminfo = h5info(stimfilename);

imgSizeX = stiminfo.Datasets.Dataspace.Size(2);%1360;
imgSizeY = stiminfo.Datasets.Dataspace.Size(3);%714;
numImgs = stiminfo.Datasets.Dataspace.Size(4);%220

nsdfolder = '/bwdata/NSDData/nsddata/experiments/nsd/';
nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
nsdDesign = load(nsdDesignFilename);
allImgs = nsdDesign.sharedix; %indices of the shared 1000 images

vecLDfolder = '/bwdata/NSDData/stimuli/vecLD';

for isub=1:8
    
    allImgs = nsdDesign.subjectim(isub,nsdDesign.masterordering);%indices of all 10000 images used for this subject
    allImgs = unique(allImgs);
    %%
    backgroundColor(1,1,:) = uint8([127,127,127]);
    
    nImgs = 1;
    
    fixPoint(1,1,:) = [255 0 0];
    
    iimg=0;
    for imgNum=allImgs
        iimg = iimg+1
        
        filename = [whichtype, 'Img' num2str(imgNum) '.mat'];
        if ~isfile(fullfile(savefolder, filename))%if file exists already no need to remake it
            imgName = ['img' num2str(imgNum) '.mat'];
            if imgNum <= 14600*1
                imgFolder = [vecLDfolder, '/images01/'];
            elseif imgNum <= 14600*2
                imgFolder = [vecLDfolder, '/images02/'];
            elseif imgNum <= 14600*3
                imgFolder = [vecLDfolder, '/images03/'];
            elseif imgNum <= 14600*4
                imgFolder = [vecLDfolder, '/images04/'];
            elseif imgNum <= 14600*5
                imgFolder = [vecLDfolder, '/images05/'];
            end

            load(fullfile(imgFolder, imgName));
           
            %% pass image through filter
            for m = 1: length(methods)
                featureMap = generateSymmetryMap(vecLD,symmetryType, methods{m}, NaN, backgroundSize, renderSize);
                model.(methods{m}) = featureMap;
            % figure;drawMATproperty(vecLD,'parallelism')
            % figure;
            % h=imagesc(featureMap);
            % set(h, 'AlphaData', 1-isnan(featureMap));
            % colormap(gca,jet);

            end

            save(fullfile(savefolder, filename),...
                'numFeatures','bandwidth','dims','model','numLevels');
        end
    end
end

