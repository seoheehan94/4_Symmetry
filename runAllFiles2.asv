%% prfSampleModel
methods = {'contour', 'medialAxis', 'area'};

% addpath(genpath('/home/hanseohe/Documents/GitHub/3rdParty'))

  % curmethod = 3;
  %   for sub = 1:8
  %       for roi = 5:7
  %           fprintf('%s. %d. %d. %s ...\n','prfSampleModel',sub,roi, methods{curmethod});
  %           prfSampleModel_symmetry(sub,roi, methods{curmethod});
  %       end
  %   end

%% regressPrfSplit

% for sub = 1:8
%     for roi =1:7
%     fprintf('%s. %d. %d. %s ...\n','regressPrfSplit',sub,roi,'contour');
%     regressPrfSplit_symmetry(sub, roi, 'par', 'contour');
%     end
% end

for sub = 1:8
    for roi =5:7
    fprintf('%s. %d. %d. %s ...\n','regressPrfSplit',sub,roi,'area');
    regressPrfSplit_symmetry(sub, roi, 'par', 'area');
    end
end

% for sub = 1:8
%     for roi =1:7
%     fprintf('%s. %d.  %s ...\n','regressPrfSplit',sub,'contour');
%     regressPrfSplit_symmetry(sub, roi, 'mir', 'contour');
%     end
% end




%% getVoxPref
% for sub = 1:8
%     fprintf('%s. %d. %d ...\n','getVoxPref_regress',sub, 1);
%     getVoxPref_regress(sub,4, 1);
% end