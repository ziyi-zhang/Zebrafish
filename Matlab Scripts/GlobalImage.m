function [res] = GlobalImage(testMat)
% Apply grid search on global 'testMat'
    
    options = optimoptions(@fminunc,'Display','none','Algorithm','quasi-newton',...
                            'MaxIterations',10,'MaxFunctionEvaluations',800);

    %% Pre-process
    testMat = double(testMat);
    testMat = testMat ./ max(testMat, [], 'all');  % normalize to 0-1   
    %% Search area
    xRange = 2:3:size(testMat, 2);
    yRange = 2:3:size(testMat, 1);
    % zRange = 16:1:24;
    zRange = 19:1:38; % for pt1
    [xArray, yArray, zArray] = meshgrid(xRange, yRange, zRange);
    xArray = xArray(:);
    yArray = yArray(:);
    zArray = zArray(:);
    fprintf('Input matrix size %d*%d*%d\n', size(testMat, 1), size(testMat, 2), size(testMat, 3));
    fprintf('X (col) search range %d to %d\n', xRange(1), xRange(end));
    fprintf('Y (row) search range %d to %d\n', yRange(1), yRange(end));
    fprintf('Z (z) search range %d to %d\n', zRange(1), zRange(end));
    %% Results
    ttlNum = length(xRange) * length(yRange) * length(zRange);
    xHist = zeros(ttlNum, 3);
    fvalHist = zeros(ttlNum, 1);
    flagHist = zeros(ttlNum, 1);
    startHist = zeros(ttlNum, 3);
    testMat = double(testMat);
    %% Search
    tic
    for count = 1:ttlNum
        x = xArray(count);
        y = yArray(count);
        z = zArray(count);
        % sub-image
        halfLength = 50;
        xmin = max(1, x-halfLength); xmax = min(size(testMat, 2), x+halfLength);
        ymin = max(1, y-halfLength); ymax = min(size(testMat, 1), y+halfLength);
        tempMat = testMat(ymin:ymax, xmin:xmax, :);
        fsigma = @(arr)GaussianSigmaFun(tempMat, arr, z);
        % quasi-newton
        [xRes, fval, exitflag, ~] = fminunc(fsigma, [x-xmin+1, y-ymin+1, 4], options);
        % store result
        xRes(1) = xRes(1) + xmin - 1;
        xRes(2) = xRes(2) + ymin - 1;
        xHist(count, :) = xRes;
        fvalHist(count) = fval;
        flagHist(count) = exitflag;
        startHist(count, :) = [x, y, z];

        fprintf('%d/%d: fval = %6.3f\n', count, ttlNum, fval);
    end
    toc

    %% make debugHist a struct
    res.x = xHist(:, 1);
    res.y = xHist(:, 2);
    res.z = startHist(:, 3);
    res.r = xHist(:, 3);
    res.fval = fvalHist;
    res.flag = flagHist;
    res.startingX = startHist(1, :);
    res.startingY = startHist(2, :);
end
