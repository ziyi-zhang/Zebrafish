function [xHist] = RefineLocation(image, points)
% Refind the location corrected by 'optical flow'
% 'image' should be the current frame
% 'points' is N-by-3 matrix

    image = double(image);
    thred = quantile(image(:), 0.97);
    image(image > thred) = thred;
    image = image ./ max(image, [], 'all');  % normalize to 0-1   

    ttlNum = size(points, 1);
    xHist = zeros(ttlNum, 3);
    fvalHist = zeros(ttlNum, 1);
    flagHist = zeros(ttlNum, 1);
    startHist = zeros(ttlNum, 3);
    % code revised from 'GlobalImage.m'
    options = optimoptions(@fminunc,'Display','none','Algorithm','quasi-newton',...
                            'MaxIterations',20,'MaxFunctionEvaluations',2000);
    for count = 1:ttlNum
        x = points(count, 1);
        y = points(count, 2);
        z = points(count, 3);
        % sub-image
        halfLength = 50;
        xmin = max(1, round(x-halfLength)); xmax = min(size(image, 2), round(x+halfLength));
        ymin = max(1, round(y-halfLength)); ymax = min(size(image, 1), round(y+halfLength));
        tempMat = image(ymin:ymax, xmin:xmax, :);
        fsigma = @(arr)GaussianSigmaFun(tempMat, arr, z);
        % quasi-newton
        [xRes, fval, exitflag, ~] = fminunc(fsigma, [x-xmin+1, y-ymin+1, 4], options);
        % store result
        xRes(1) = xRes(1) + xmin - 1;
        xRes(2) = xRes(2) + ymin - 1;
        xHist(count, :) = [xRes(1:2), z];  % xRes(3) is radius!
        fvalHist(count) = fval;
        flagHist(count) = exitflag;
        startHist(count, :) = [x, y, z];

        fprintf('%d/%d: fval = %6.3f\n', count, ttlNum, fval);
    end

end
