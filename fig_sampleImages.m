%%
clear all;

img_file = '/bwlab/Users/SeoheeHan/NSDData/rothzn/nsd/stimuli/mirfilter/mirImg28.mat';
load(img_file);
mainFileName = 'img28';

%img size
topLeft = [77, 77];
imgWidth =360;
imgHeight = 360;

model.area = imcrop(model.area, [topLeft(1), topLeft(2), imgHeight-1, imgWidth-1]);
model.contour = imcrop(model.contour, [topLeft(1), topLeft(2), imgHeight-1, imgWidth-1]);
model.medialAxis = imcrop(model.medialAxis, [topLeft(1), topLeft(2), imgHeight-1, imgWidth-1]);


%colormap
base_cmap = turbo(numCols);  
cmap = [1 1 1;base_cmap];

% area
fig = figure; 
imagesc(model.area);
set(gca, 'XTick', [], 'YTick', []);
box on;colormap(cmap);
axis image;
frame = getframe(gca);
img_area = frame.cdata;
img_area = imresize(img_area, [imgHeight, imgWidth]);
% figure; imshow(LD_img)
imwrite(img_area, [mainFileName, '_area.png']);

% contour
fig = figure; 
imagesc(model.contour);
set(gca, 'XTick', [], 'YTick', []);
box on;colormap(cmap);
axis image;
frame = getframe(gca);
img_contour = frame.cdata;
img_contour = imresize(img_contour, [imgHeight, imgWidth]);
% figure; imshow(LD_img)
imwrite(img_contour, [mainFileName, '_contour.png']);

% medialAxis
fig = figure; 
imagesc(model.medialAxis);
set(gca, 'XTick', [], 'YTick', []);
box on;colormap(cmap);
axis image;
frame = getframe(gca);
img_medialAxis = frame.cdata;
img_medialAxis = imresize(img_medialAxis, [imgHeight, imgWidth]);
% figure; imshow(LD_img)
imwrite(img_medialAxis, [mainFileName, '_medialAxis.png']);

  

