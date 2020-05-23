function [image_gau] = GaussianFilter(image, hsize, sigma, showImage)
% [image_gau] = GaussianFilter(image, hsize, sigma, showImage)
% Apply Gaussian filter on image.
% If 'image' is 3D array, it will be processed as a stack of 2D images
% Ziyi. Feb, 2020.

    if nargin<4, showImage=false;end
    if nargin<3 || isempty(sigma), sigma = 3;end
    if nargin<2 || isempty(hsize), hsize = 8;end

    h = fspecial('gaussian', hsize, sigma);
    image_gau = zeros(size(image));
    for i = 1:size(image, 3)
        image_gau(:, :, i) = imfilter(image(:, :, i), h);
    end

    if showImage
        figure
        image = image_gau(:, :, 1);
        imshow(image, [min(image, [], 'all'), max(image, [], 'all')]);
        colormap bone
        titleStr = "l:" + hsize + " s:" + sigma;
        title(titleStr);
    end
end
