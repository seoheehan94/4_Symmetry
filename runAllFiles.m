%% prfSampleModel
methods = {'contour', 'medialAxis', 'area'};

% addpath(genpath('/home/hanseohe/Documents/GitHub/stimulusVignetting'))
% curmethod = 1;
% % for curmethod = 1:3
%     for sub = 7:8
%         for roi = 5:7
%             fprintf('%s. %d. %d. %s ...\n','prfSampleModel',sub,roi, methods{curmethod});
%             prfSampleModel_symmetry(sub,roi, methods{curmethod});
%         end
%     end
% % end

% curmethod = 2;
%     % for sub = 4:8
%     %     for roi = 2
%     %         fprintf('%s. %d. %d. %s ...\n','prfSampleModel',sub,roi, methods{curmethod});
%     %         prfSampleModel_symmetry(sub,roi, methods{curmethod});
%     %     end
%     % end
% 
%     for sub = 6:8
%         for roi = 5:7
%             fprintf('%s. %d. %d. %s ...\n','prfSampleModel',sub,roi, methods{curmethod});
%             prfSampleModel_symmetry(sub,roi, methods{curmethod});
%         end
%     end


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
for sub = 1:8
    fprintf('%s. %d. %d. %s ...\n','getVoxPref_symmetry',sub, 7, 'area');
    getVoxPref_symmetry(sub,7,'par', 'area')
end