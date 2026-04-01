clear all;
% methods = {'contour', 'medialAxis', 'area'};
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
% whichmethod='contour';
prompt = "Type of method?";
whichmethod = input(prompt,"s");

filefolder = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/';
savefolder = ['/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/',whichmethod,'Allfilter/'];%to save model outputs

numLevels = 1;
numFeatures = 3;
bandwidth=1;
dims=[512,512];

for imgNum = 1:73000
    clearvars model filename par mir tap parfile mirfile tapfile
    filename = [whichmethod, 'Img' num2str(imgNum) '.mat'];
    fprintf('%s ...\n',filename);
    parfile = ['parfilter/parImg' num2str(imgNum) '.mat'];
    mirfile = ['mirfilter/mirImg' num2str(imgNum) '.mat'];
    tapfile = ['tapfilter/tapImg' num2str(imgNum) '.mat'];
    par=load([filefolder, parfile]);
    mir=load([filefolder, mirfile]);
    tap=load([filefolder, tapfile]);

    % model(3x512x512)
    
    model(1,:,:)=par.model.(whichmethod);
    model(2,:,:)=mir.model.(whichmethod);
    model(3,:,:)=tap.model.(whichmethod);


    save(fullfile(savefolder, filename),...
        'numFeatures','bandwidth','dims','model','numLevels');

end