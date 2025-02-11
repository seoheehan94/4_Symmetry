%% prfSampleModel
methods = {'contour', 'medialAxis', 'area'};

addpath(genpath('/home/hanseohe/Documents/GitHub/3rdParty'))
    for sub = 5:8
            fprintf('%s. %d. %d. %s ...\n','prfSampleModel',sub,1, methods{2});
            prfSampleModel_symmetry(sub,1, methods{2});
    end
    for sub = 1:8
            fprintf('%s. %d. %d. %s ...\n','prfSampleModel',sub,1, methods{3});
            prfSampleModel_symmetry(sub,1, methods{3});
    end

%% regressPrfSplit
% for sub = 1:8
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 1);
% regressPrfSplit(sub, 1);
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 2);
% regressPrfSplit(sub, 2);
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 3);
% regressPrfSplit(sub, 3);
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 4);
% regressPrfSplit(sub, 4);
% end
% % 
% for sub = 1:8
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 1);
% regressPrfSplit_control(sub, 1);
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 2);
% regressPrfSplit_control(sub, 2);
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 3);
% regressPrfSplit_control(sub, 3);
% fprintf('%s. %d. %d ...\n','regressPrfSplit',sub, 4);
% regressPrfSplit_control(sub, 4);
% end

%% getVoxPref
% for sub = 1:8
%     fprintf('%s. %d. %d ...\n','getVoxPref_regress',sub, 1);
%     getVoxPref_regress(sub,4, 1);
% end