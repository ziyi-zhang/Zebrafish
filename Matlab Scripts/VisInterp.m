function [] = VisInterp(mat, dim, l, supersample)
%

    if nargin<4
        supersample = 4;
    end
    mat = double(mat);

    [sx, sy, sz] = size(mat);
    if dim==1

    elseif dim==2

    elseif dim==3
        res = zeros(sx*supersample, sy*supersample);
        [xArray, yArray] = meshgrid(1/supersample:1/supersample:sx, 1/supersample:1/supersample:sy);
        xArray_ = xArray(:);
        yArray_ = yArray(:);
        lArray = repmat(double(l), length(xArray_), 1);
        resArray = interp3(mat, xArray_, yArray_, lArray, 'spline');
        for i = 1:size(xArray, 2)
            for j = 1:size(xArray, 1)
                res(j, i) = resArray((i-1)*size(xArray, 1)+j);
            end
        end
    end

    figure
    subplot(1, 2, 1);
    if dim==1
        mat = mat(l, :, :);
        mat = squeeze(mat);
    elseif dim==2
        mat = mat(:, l, :);
        mat = squeeze(mat);
    elseif dim==3
        mat = mat(:, :, l);
        mat = squeeze(mat);
    end
    imshow(mat, [min(mat, [], 'all'), max(mat, [], 'all')]);
    colormap('jet');
    colorbar
    titleStr = sprintf("Original Image z=%d", l);
    title(titleStr)
    subplot(1, 2, 2);
    imshow(res, [quantile(res(:), 0.05), quantile(res(:), 0.95)]);
    %imshow(res, [0.7e4, 2.4e4]);
    %imshow(res, [min(res, [], 'all'), max(res, [], 'all')])
    %imagesc(res)
    colormap('jet');
    colorbar
    titleStr = sprintf("Interpolated image with supersample=%dx", supersample);
    title(titleStr)
end
