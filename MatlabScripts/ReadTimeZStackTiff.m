function [image_] = ReadTimeZStackTiff(filename)
% [image_] = ReadTimeZStackTiff(filename)
% read a TIFF file and return a struct with subimages, corresponding
% metadata and the TIFF imageInfo
% Ziyi. Nov, 2020.
    
    if ~strcmp(filename(end), 'f')
        filename = strcat(filename, '.tif');
    end
    % Read a multipage tiff, assuming each file is the same size
    t = Tiff(filename, 'r');
    image(:,:,1) = t.read(); % Read the first image to get the array dimensions correct.
    if t.lastDirectory()
         image_ = image;
         return; % If the file only contains one page, we do not need to continue.
    end
    % Read all remaining pages (directories) in the file
    t.nextDirectory();
    while true
        image(:, :, end+1) = t.read(); %#ok
        if t.lastDirectory()
            break;
        else
            t.nextDirectory();
        end
    end
    imageInfo = imfinfo(filename);
    
    % now label images
    description = imageInfo(1).ImageDescription;
    totalImages = extractNumber(description, 'images=');
    slices = extractNumber(description, 'slices=');
    frames = extractNumber(description, 'frames=');

    meta = zeros(2, totalImages); % channel, slice, frame
    count = 0;
    for i = 1:frames
        for j = 1:slices
            count = count + 1;
            meta(:, count) = [j; i];
        end
    end
    
    % composite
    image_.image = image;
    image_.meta = meta;
    image_.info = imageInfo;
    
    % print summary
    fprintf('The TIFF image <%s> has %d subimages.\n', filename, totalImages);
    fprintf('Slices=%d, Frames=%d\n', slices, frames);
    fprintf('Each subimage has width=%d, height=%d, bitDepth=%d\n', imageInfo(1).Width, imageInfo(1).Height, imageInfo(1).BitDepth);
end


function [res] = extractNumber(charArr, key)
% return the number corresponds to the 'key'
% -1 on exception

    loc = strfind(charArr, key);
    if (isempty(loc) || length(loc) > 1)
        res = -1;
        return;
    end
    start_ = loc + length(key);
    end_ = start_;
    while (end_ <= length(charArr) && charArr(end_)>='0' && charArr(end_)<='9')
        end_ = end_ + 1;
    end
    if (end_-start_ == 0)
        res = -1;
        return;
    end
    res = str2double(charArr(start_:end_-1));
end
