function [image_log] = LoGFilter(image, hsize, sigma, showImage)
% [image_log] = LoGFilter(image, hsize, sigma, showImage)
% Apply LoG (Laplacian of Gaussian) filter on image.
% If 'image' is 3D array, it will be processed as a stack of 2D images
% Ziyi. Feb, 2020.

    if nargin<4, showImage=false;end
    if nargin<3 || isempty(sigma), sigma = 2;end
    if nargin<2 || isempty(hsize), hsize = 12;end

    h = fspecial('log', hsize, sigma);
    image_log = zeros(size(image));
    for i = 1:size(image, 3)
        image_log(:, :, i) = imfilter(image(:, :, i), h);
    end

    if showImage
        figure
        image = image_log(:, :, 1);
        imshow(image, [min(image, [], 'all'), max(image, [], 'all')]);
        colormap bone
        titleStr = "l:" + hsize + " s:" + sigma;
        title(titleStr);
    end
end
