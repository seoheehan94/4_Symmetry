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

    curpar=parModel.model.contour;
    curmir=mirModel.model.contour;

        if length(nonzeros(curpar)) == length(nonzeros(curmir))
            curpar_nonzeros=nonzeros(curpar);
            curmir_nonzeros=nonzeros(curmir);
            
        else
            mask = (curpar == 0 & curmir ~= 0) | (curmir == 0 & curpar ~= 0);
            % Replace all 0s in both A and B with NaN
            curpar(curpar == 0) = NaN;
            curmir(curmir == 0) = NaN;
            % Restore 0s where one was originally 0 but the other was not
            curpar(mask) = 0;
            curmir(mask) = 0;
            curpar_nonzeros=curpar(~isnan(curpar));
            curmir_nonzeros=curmir(~isnan(curmir));
        end

        [R,P] = corrcoef(curpar_nonzeros,curmir_nonzeros);
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