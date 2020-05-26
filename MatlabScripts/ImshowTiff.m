function [] = ImshowTiff(image, index, enhanced, zoomIn)
% [] = ImshowTiff(image[, index][, enhanced][, zoomIn])
% [] = ImshowTiff(mat)
% 'image' should be a struct with meta data created by 'readMultipageTiff.m'
% 'image.image' is a 3D matrix of size [length*width*stacks] eg:
% [1024*1024*1558 uint16]
% 'index' is an array containing the desired indexes of subimages
% 'enhanced' is a flag indicating whether the output should be shown with 
% enhanced contrast
% 'zoomIn' is the ratio of zoom in and appies on all subplots
% 'mat' is a 2D/3D matrix. The contrast will be auto enhanced. If the input
% is 3D, the third dimension will be summed
% Ziyi. Feb, 2020.

    figure
    if size(image, 1)>1 % if 2D/3D matrix
        if size(image, 3) > 1
            image = sum(image, 3);
        end
        imshow(image, [min(image, [], 'all'), max(image, [], 'all')]);
        colormap(bone); % use bone instead of pink
        return;
    end
    if ~isstruct(image)
        warning("Error using imshowTiff.");
        return;
    end
    
    if nargin<4, zoomIn = 1;end
    if nargin<3, enhanced = false;end
    if nargin<2 || isempty(index), index = 1:size(image.image, 3);end
    N = length(index);
    ncol = ceil(sqrt(16/9*N));
    nrow = ceil(N / ncol);
    if (N==1) ncol=1;nrow=1;end
    % zoom in calculation
    imageX = size(image.image, 1);
    imageY = size(image.image, 2);
    xzoomIn(1) = ceil(imageX / 2 * (1 - 1.0/zoomIn));
    xzoomIn(2) = imageX - xzoomIn(1);
    yzoomIn(1) = ceil(imageY / 2 * (1 - 1.0/zoomIn));
    yzoomIn(2) = imageY - yzoomIn(1);
    if xzoomIn(1)<1, xzoomIn(1)=1;end
    if yzoomIn(1)<1, yzoomIn(1)=1;end
    % enhance contrast
    if enhanced
        t = image.image(:, :, index);
        t = t(:);
        lowerBound = quantile(t, 0.003);
        upperBound = quantile(t, 0.997);
    end
    % plot
    ha = tight_subplot(nrow, ncol);
    for i = 1:N

        axes(ha(i)); %#ok
        if (enhanced)
            imshow(image.image(xzoomIn(1):xzoomIn(2), yzoomIn(1):yzoomIn(2), index(i)), [lowerBound, upperBound]);
        else
            imshow(image.image(xzoomIn(1):xzoomIn(2), yzoomIn(1):yzoomIn(2), index(i)));
        end
        colormap(pink);
        if (isfield(image, 'meta'))
            titleStr = "C:" + image.meta(1, index(i)) + " S:" + image.meta(2, index(i)) + " F:" + image.meta(3, index(i));
            title(titleStr);
        else 
            title(i);
        end 
    end
end