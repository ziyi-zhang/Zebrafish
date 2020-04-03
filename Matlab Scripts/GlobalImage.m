function [xHist, fvalHist, flagHist, debugHist] = GlobalImage(testMat)
% Apply grid search on global 'testMat'
    
    options = optimoptions(@fminunc,'Display','none','Algorithm','quasi-newton',...
                            'MaxIterations',20,'MaxFunctionEvaluations',800);
                        
    %% Search area
    xRange = 5:3:size(testMat, 2);
    yRange = 5:3:size(testMat, 1);
    % zRange = 16:3:24;
    zRange = 19:3:38;
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
    xHist = zeros(ttlNum, 4);
    fvalHist = zeros(ttlNum, 1);
    flagHist = zeros(ttlNum, 1);
    debugHist = zeros(ttlNum, 3);
    testMat = double(testMat);
    %% Search
    tic
    parfor count = 1:ttlNum
        x = xArray(count);
        y = yArray(count);
        z = zArray(count);
        % sub-image
        halfLength = 15;
        xmin = max(1, x-halfLength); xmax = min(size(testMat, 2), x+halfLength);
        ymin = max(1, y-halfLength); ymax = min(size(testMat, 1), y+halfLength);
        tempMat = testMat(ymin:ymax, xmin:xmax, :);
        fsigma = @(arr)GaussianSigmaFun(tempMat, arr);
        % quasi-newton
        [res, fval, exitflag, ~] = fminunc(fsigma, [x-xmin+1, y-ymin+1, z, 4], options);
        % store result
        res(1) = res(1) + xmin - 1;
        res(2) = res(2) + ymin - 1;
        xHist(count, :) = res;
        fvalHist(count) = fval;
        flagHist(count) = exitflag;
        debugHist(count, :) = [x, y, z];

        fprintf('%d/%d: fval = %6.3f\n', count, ttlNum, fval);
    end
    toc
end
