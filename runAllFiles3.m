%% prfSampleModel
methods = {'contour', 'medialAxis', 'area'};
% 
% addpath(genpath('/home/hanseohe/Documents/GitHub/3rdParty'))
% 
% curmethod = 3;
%     for sub = 1:8
%         for roi = 2:4
%             fprintf('%s. %d. %d. %s ...\n','prfSampleModel',sub,roi, methods{curmethod});
%             prfSampleModel_symmetry(sub,roi, methods{curmethod});
%         end
%     end


    % for sub = 1:8
    %         fprintf('%s. %d. %d. %s ...\n','prfSampleModel',sub,1, methods{3});
    %         prfSampleModel_symmetry(sub,1, methods{3});
    % end
    % 
%% regressPrfSplit
for sub = 1:8
fprintf('%s. %d.  %s ...\n','regressPrfSplit',sub, 'medialAxis');
regressPrfSplit_symmetry(sub, [1,2,3,4,5,6,7], 'par', 'medialAxis');
end

%% getVoxPref
% for sub = 1:8
%     fprintf('%s. %d. %d ...\n','getVoxPref_symmetry',sub, 7);
%     getVoxPref_symmetry(sub,7,'medialAxis')
% end
% for sub = 1:8
%     fprintf('%s. %d. %d ...\n','getVoxPref_symmetry',sub, 7);
%     getVoxPref_symmetry(sub,7,'area')
% end
