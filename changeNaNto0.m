methods = {'contour', 'medialAxis', 'area'};
% cd '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli';
%pyramidfolder = '/misc/data18/rothzn/nsd/stimuli/pyramid/';%to save model outputs
% switch (lower(symmetryType))
%     case 'parallelism'
%         whichtype = 'par';
%     case 'separation'
%         whichtype = 'sep';
%     case 'mirror'
%         whichtype = 'mir';
%     case 'taper'
%         whichtype = 'tap';
% end
whichtype='tap';
savefolder = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/',whichtype,'filter/'];%to save model outputs

for imgNum = 1:73000
    filename = [whichtype, 'Img' num2str(imgNum) '.mat'];
    fprintf('%s ...\n',filename);
    load([savefolder, filename]);
    for curMethod = 1:3
        model.(methods{curMethod})(isnan(model.(methods{curMethod}))) = 0;


        save(fullfile(savefolder, filename),...
            'numFeatures','bandwidth','dims','model','numLevels');

    end

end