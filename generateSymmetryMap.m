function featureMap = generateSymmetryMap(vecLD, symmetryType, method, backgroundValue, backgroundSize,renderSize)
%   The list of options include:
%   1. 'parallelism'
%   2. 'separation'
%   3. 'mirror'
%   4. 'taper'
if nargin < 6
    renderSize = vecLD.imsize;
end
if nargin < 5
    backgroundSize = renderSize;
end
if nargin < 4
    backgroundValue = NaN;
end
featureMap = zeros(vecLD.imsize) + backgroundValue;

switch (lower(symmetryType))
    case 'parallelism'
        allScores = vecLD.parallelism_allScores;
        allY = vecLD.parallelism_allY;
        allX = vecLD.parallelism_allX;
    case 'separation'
        allScores = vecLD.separation_allScores;
        allY = vecLD.separation_allY;
        allX = vecLD.separation_allX;
    case 'mirror'
        allScores = vecLD.mirror_allScores;
        allY = vecLD.mirror_allY;
        allX = vecLD.mirror_allX;
    case 'taper'
        allScores = vecLD.taper_allScores;
        allY = vecLD.taper_allY;
        allX = vecLD.taper_allX;
    otherwise
        error(['Unknown symmetryType string: ',symmetryType]);
end

switch (lower(method))
    case 'contour'

        for p = 1:numel(allScores)
            featureMap(allY(p),allX(p)) = allScores(p);
        end

        % % Resize the feature map to the render size
        % resizedFeatureMap = imresize(featureMap, renderSize, 'nearest');
        %

    case 'medialaxis'
        img = renderLinedrawing(vecLD);
        MAT = computeMAT(img);
        skeletalBranches = traceSkeleton(MAT);
        [skeletonImageWithRating,skeletalBranches] = computeMATproperty(MAT,symmetryType,skeletalBranches);
        featureMap = skeletonImageWithRating;
        featureMap(featureMap==0)=NaN;
    case 'area'
        img = renderLinedrawing(vecLD);
        MAT = computeMAT(img);
        skeletalBranches = traceSkeleton(MAT);
        [skeletonImageWithRating,skeletalBranches] = computeMATproperty(MAT,symmetryType,skeletalBranches);
        contourImageWithRating = mapMATtoContourFullTangential(skeletalBranches, img, skeletonImageWithRating);
        contourImageWithRating(skeletonImageWithRating~=0)=skeletonImageWithRating(skeletonImageWithRating~=0);
        contourImageWithRating(contourImageWithRating==0)=NaN;
        featureMap = contourImageWithRating;

    otherwise
        error(['Unknown method string: ',method]);
end


resizedFeatureMap = maxPoolResize(featureMap, renderSize);

scaleY = renderSize(1)/vecLD.imsize(1);
scaleX = renderSize(2)/vecLD.imsize(2);
offsetY = round((backgroundSize(1) - renderSize(1))/2);
offsetX = round((backgroundSize(2) - renderSize(2))/2);


% % Initialize a padded feature map with the background value
paddedFeatureMap = zeros(backgroundSize) + backgroundValue;

% Place the resized map in the center of the padded map
paddedFeatureMap(offsetY + 1:offsetY + renderSize(1), offsetX + 1:offsetX + renderSize(2)) = resizedFeatureMap;

% Return the padded feature map
featureMap = paddedFeatureMap;

    function resizedFeatureMap = maxPoolResize(featureMap, renderSize)
        [origHeight, origWidth] = size(featureMap);
        newHeight= renderSize(1);
        newWidth= renderSize(2);

        scaleY = origHeight / newHeight;
        scaleX = origWidth / newWidth;

        resizedFeatureMap = zeros(newHeight, newWidth) * NaN;

        for y = 1:newHeight
            for x = 1:newWidth
                % Calculate the original region corresponding to this pixel
                yStart = floor((y - 1) * scaleY) + 1;
                yEnd = min(floor(y * scaleY), origHeight);
                xStart = floor((x - 1) * scaleX) + 1;
                xEnd = min(floor(x * scaleX), origWidth);

                % Take the maximum value in the region
                region = featureMap(yStart:yEnd, xStart:xEnd);
                resizedFeatureMap(y, x) = max(region(:), [], 'omitnan');
            end
        end
    end



end