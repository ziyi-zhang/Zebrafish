function [image, imageInfo] = ReadZstackTiff(filename)
% [image, imageInfo] = ReadZstackTiff(filename)
% read a TIFF file and return a struct with subimages
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
end
