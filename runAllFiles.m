%% regressPrfSplit
sym_type = {'mir', 'tap', 'par'};
% % 
% for sub = 1:8
%     for roi =1:7
%     fprintf('%s. %d. %d. %s ...\n','regressPrfSplit',sub,roi,'par');
%     regressPrfSplit_symmetry_new(sub, roi, 'par');
%     end
% end

for sub = 1:8
    for roi = 8 
    fprintf('%s. %d. %d. %s ...\n','regressPrfSplit',sub,roi,'par');
    regressPrfSplit_symmetry_new(sub, roi, 'par');
    end
end


for sub = 1:8

    getVoxPref_symmetry_new(sub, 1:8, 'par')
end