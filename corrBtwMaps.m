close all;
clear all;


parfolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/parfilter/';%to save model outputs
parList = dir(parfolder);
parList = parList(~ismember({parList(:).name},{'.','..'}));

mirfolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/parfilter/';%to save model outputs
mirList = dir(mirfolder);
mirList = mirList(~ismember({mirList(:).name},{'.','..'}));


savefolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/corrBtwMaps/';%to save model outputs

%%

iimg=0;
for imgNum=1:length(parList)
    iimg = iimg+1
    parname = parList(imgNum).name;
    mirname = mirList(imgNum).name;

    filename = ['maImg' num2str(imgNum) '.mat'];
    if ~isfile(fullfile(savefolder, filename))%if file exists already no need to remake it

        parModel = load(fullfile(parfolder, parname), 'model');
        mirModel = load(fullfile(mirfolder, mirname), 'model');

        [R,P] = corrcoef(parModel.model.contour,mirModel.model.contour);


        
        save(fullfile(savefolder, filename),...
            'numFeatures','bandwidth','dims','model','numLevels');
    end
end

