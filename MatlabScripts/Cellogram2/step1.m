% data read in
% [image, cropData] = step1('C:\Users\ziyiz\Downloads\Data-TFM\Data-TFM\20.11.14-9', 5, 11, 871, 574, 1096, 760);

% [image, cropData] = step1('C:\Users\ziyiz\Downloads\Data-TFM\Data-TFM\20.10.23-8\20.10.23-8-1-for', 2, 27, 353, 605, 749, 939);
% [image, cropData] = step1('C:\Users\ziyiz\Downloads\Data-TFM\Data-TFM\20.10.17-8\20.10.17-8-1-for', 1, 25, 852, 924, 1138, 1200);

function [image_, cropData] = step1(fileName, sliceBegin, sliceEnd, r0, c0, r1, c1)

    image = ReadTimeZStackTiff(fileName);
    
    % want all frames, slices in [sliceBegin, sliceEnd]
    if (~isempty(image.meta))
        rawSlices = max(image.meta(1, :));
        frames = max(image.meta(2, :));
    else
        rawSlices = size(image.image, 3);
        frames = 1;
    end
    image_ = cell(frames, 1);
    
    for f = 1:frames
        baseIdx = (f-1) * rawSlices;
        % normalize each frame independently
        img = double(image.image(r0:r1, c0:c1, baseIdx+sliceBegin:baseIdx+sliceEnd));
        minPixel = min(img, [], 'all');
        % maxPixel = max(img, [], 'all');
        maxPixel = quantile(img(:), 0.992);
        img = (img - minPixel) ./ (maxPixel - minPixel);
        img(img > 1) = 1;
        image_{f} = img;
    end

    cropData.sliceBegin = sliceBegin;
    cropData.sliceEnd = sliceEnd;
    cropData.r0 = r0;
    cropData.c0 = c0;
    cropData.r1 = r1;
    cropData.c1 = c1;
end
