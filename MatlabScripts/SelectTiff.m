function [image] = SelectTiff(image, channel, slice, frame, indexStart, indexEnd)
% [image] = SelectTiff(image[, channel][, slice][, frame][, indexStart][, indexEnd])
% select a subset of subimages from 'image' with five constraints
% 'channel', 'slice', 'frame', 'indexStart' and 'indexEnd'
% Example Use:
% selectTiff(image, 1, [], 1) % select all subimages from first frame with
% channel one
% selectTiff(image, 1, 1:5:41, [], 1, 300) % select from first 300
% subimages with a fixed interval 5
% Ziyi. Feb, 2020.

    % process input
    if nargin < 6 || isempty(indexEnd)
        indexEnd = size(image.image, 3);
    end
    if nargin < 5 || isempty(indexStart)
        indexStart = 1;
    end
    if nargin < 4
        frame = [];
    end
    if nargin < 3
        slice = [];
    end
    if nargin < 2
        channel = [];
    end
    
    % select
    keepMask = false(1, size(image.image, 3));
    for i = indexStart:indexEnd
        if ~isempty(channel)
            if ~ismember(image.meta(1, i), channel)
                continue;
            end
        end
        if ~isempty(slice)
            if ~ismember(image.meta(2, i), slice)
                continue;
            end
        end
        if ~isempty(frame)
            if ~ismember(image.meta(3, i), frame)
                continue;
            end
        end
        keepMask(i) = true;
    end
    image.image = image.image(:, :, keepMask);
    image.meta = image.meta(:, keepMask);
    
    % log
    fprintf('Selected %d subimages from %d subimages.\n', nnz(keepMask), length(keepMask));
end
