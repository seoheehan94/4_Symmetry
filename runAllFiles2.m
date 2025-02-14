%% prfSampleModel
methods = {'contour', 'medialAxis', 'area'};

% addpath(genpath('/home/hanseohe/Documents/GitHub/3rdParty'))

  curmethod = 3;
    for sub = 1:8
        for roi = 5:7
            fprintf('%s. %d. %d. %s ...\n','prfSampleModel',sub,roi, methods{curmethod});
            prfSampleModel_symmetry(sub,roi, methods{curmethod});
        end
    end

%% regressPrfSplit
for sub = 1:8
fprintf('%s. %d.  %s ...\n','regressPrfSplit',sub, 'contour');
regressPrfSplit_symmetry(sub, [1,2,3,4,5,6], 'contour');

end
% 
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