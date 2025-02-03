% compute Taper batch run code
clear all;
rng(4228);

% Add path
addpath(genpath('/home/hanseohe/Documents/GitHub/MLV_toolbox'));
original_folder = '/bwdata/NSDData/stimuli/vecLD'; % Original image file folder
addpath(original_folder);
original_subfolder = dir(original_folder);
%original_subfolder = original_subfolder(~ismember({original_subfolder(:).name},{'.','..'}));
original_subfolder = original_subfolder(ismember({original_subfolder(:).name},{'images01','images02','images03','images04','images05'}));

save_folder = '/bwdata/NSDData/stimuli/vecLD';

%% Read files
for k = 1 : length(original_subfolder)

    name = original_subfolder(k).name;
    imFilePath = natsortfiles(dir(strcat(original_folder, '/', name,'/', '*.mat')));

    save_subfolder = fullfile(save_folder, name);


    for j = 1 : length(imFilePath)
        vecLD = [];
        imFile = imFilePath(j).name;
        fprintf('%d. %s ...\n',j,imFile);
        [~,mainFileName,~] = fileparts(imFile);
        filePath = strcat(original_folder, '/', name,'/', imFile);
        load(filePath);

        %% This is the actual process for a single vecLD
        img = renderLinedrawing(vecLD);
        MAT = computeMAT(img);
        [MATimg,MATskel,branches] = computeAllMATproperties(MAT,img, {'taper'});
        properties = fieldnames(MATimg);

        for p = 1:length(properties)
            thisPropImg = mapMATtoContour(branches,img,MATskel.(properties{p}));
            vecLD = MATpropertiesToContours(vecLD,thisPropImg,properties{p});
            vecLD = getMATpropertyStats(vecLD,properties{p});
        end
        %% Save results
        save(strcat(save_subfolder, '/', mainFileName, '.mat'), 'vecLD')

    end
end