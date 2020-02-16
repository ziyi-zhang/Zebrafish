function [image] = LoGFilter(image, hsize, sigma, showImage)
% [image] = LoGFilter(image, hsize, sigma, showImage)
% Apply LoG (Laplacian of Guassian) filter on image.

    if nargin<4, showImage=false;end
    if nargin<3 || isempty(sigma), sigma = 2;end
    if nargin<2 || isempty(hsize), hsize = 8;end

    h = fspecial('log', hsize, sigma);
    image = imfilter(image, h);
    if showImage
        figure
        imshow(image, [min(image, [], 'all'), max(image, [], 'all')]);
        colormap bone
        titleStr = "l:" + hsize + " s:" + sigma;
        title(titleStr);
    end
end
