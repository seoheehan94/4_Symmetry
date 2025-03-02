close all;
clear all;


parfolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/parfilter/';%to save model outputs
parList = dir(parfolder);
parList = parList(~ismember({parList(:).name},{'.','..'}));

mirfolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/mirfilter/';%to save model outputs
mirList = dir(mirfolder);
mirList = mirList(~ismember({mirList(:).name},{'.','..'}));


savefolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/Symmetry/corrBtwMaps/';%to save model outputs
filePath = fullfile(savefolder, 'results_corrBtwMaps_nonzeros.mat');
if isfile(filePath)
    load(filePath, 'allResultsTable');
end
%%

iimg=0;
for imgNum=1:length(parList)
    iimg = iimg+1
    parname = parList(imgNum).name;
    mirname = mirList(imgNum).name;


        parModel = load(fullfile(parfolder, parname), 'model');
        mirModel = load(fullfile(mirfolder, mirname), 'model');

        [R,P] = corrcoef(nonzeros(parModel.model.contour),nonzeros(mirModel.model.contour));

        resultsTable = table(imgNum, R(1,2), P(1,2), ...
            'VariableNames', {'imgnum', 'R', 'P'});

        if exist('allResultsTable', 'var')
            allResultsTable = [allResultsTable; resultsTable];
        else
            allResultsTable = resultsTable;
        end
    save(filePath, 'allResultsTable');
end

% writetable(allResultsTable, 'results_corrBtwMaps.csv');