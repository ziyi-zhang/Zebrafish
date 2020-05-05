function [points, score] = ScorePoints2(inStruct, visPlot, image)
% After running 'GlobalImage', score the resultant points

    if nargin<2
        visPlot = 0;
    end
    
    % filter config
    fvalThres = -0.1;
    radiusThres = 7;
    membraneThres = 0.7; % binary threshold
    HistPeriThres = 0.85; % quantile threshold (of peripheral samples' variance)
    % cluster config
    thresDist = 0.45;
    thresMerge = 10;
    thresClusterSize = 8;
    
    % retrive data
    x = inStruct.x;
    y = inStruct.y;
    z = inStruct.z;
    r = inStruct.r;
    fval = inStruct.fval;
    flag = inStruct.flag;
    startingX = inStruct.startingX;
    startingY = inStruct.startingY;
    fprintf("#starting points                = %d\n", length(x));
    
    % want fval < fvalThres
    mask = fval < fvalThres;
    Filter(mask);
    fprintf("After filter fval:      #points = %d\n", length(x));
    
    % want radius <= radiusThres
    mask = r <= radiusThres;
    Filter(mask);
    fprintf("After filter radius:    #points = %d\n", length(x));
    
    % want peri color bright (in membrane area)
    biImage = image;
    for i = 1:size(biImage, 3)
        t = biImage(:, :, i);
        t = t(:);
        biImage(:, :, i) = biImage(:, :, i) > quantile(t, membraneThres);
    end
    count = interp3(biImage, x, y+r*1.2, z, 'nearest')+...
            interp3(biImage, x, y-r*1.2, z, 'nearest')+...
            interp3(biImage, x+r*1.2, y, z, 'nearest')+...
            interp3(biImage, x-r*1.2, y, z, 'nearest')+...
            interp3(biImage, x+r*0.85, y+r*0.85, z, 'nearest')+...
            interp3(biImage, x+r*0.85, y-r*0.85, z, 'nearest')+...
            interp3(biImage, x-r*0.85, y+r*0.85, z, 'nearest')+...
            interp3(biImage, x-r*0.85, y-r*0.85, z, 'nearest');

    mask = count >= 8;
    Filter(mask);
    fprintf("After filter membrane:  #points = %d\n", length(x));
    
    % hist of peripheral area
    histPeri = [interp3(image, x, y+r*1.2, z, 'nearest'),...
                interp3(image, x, y-r*1.2, z, 'nearest'),...
                interp3(image, x+r*1.2, y, z, 'nearest'),...
                interp3(image, x-r*1.2, y, z, 'nearest'),...
                interp3(image, x+r*0.85, y+r*0.85, z, 'nearest'),...
                interp3(image, x+r*0.85, y-r*0.85, z, 'nearest'),...
                interp3(image, x-r*0.85, y+r*0.85, z, 'nearest'),...
                interp3(image, x-r*0.85, y-r*0.85, z, 'nearest')];
    histPeriVar = var(double(histPeri), 0, 2);
    mask = histPeriVar < quantile(histPeriVar, HistPeriThres);
    Filter(mask);
    fprintf("After filter histogram: #points = %d\n", length(x));
    
    % create a struct
    points = table;
    points.x = x;
    points.y = y;
    points.z = z;
    points.r = r;
    points.fval = fval;
    points.flag = flag;
    points.startingX = startingX;
    points.startingY = startingY;
    
    %% cluster
    count = 0;
    idx = zeros(1, height(points));
    maski = false(0, height(points));
    for i = 1:height(points)
    
        if (idx(i) > 0), continue;end
        xy = table2array(points(:, 1:2));
        xy = xy - xy(i, :);
        dist = xy(:, 1).^2 + xy(:, 2).^2;
        mask = dist < thresDist^2;
        if (nnz(mask) < thresClusterSize), continue;end
        % new cluster
        count = count + 1;
        maski(count, :) = mask;
        idx(mask) = count;
    end
    % count size again
    clusterSize = zeros(count, 1);
    meanXYZR = zeros(count, 4);
    for i = 1:count
        mask = idx == i;
        doubleArray = table2array(points(mask, :));
        clusterSize(i) = nnz(mask);
        meanXYZR(i, :) = mean(doubleArray(:, 1:4));
    end
    maskSize = clusterSize > thresClusterSize;
    % merge
    count = 0;
    idx2 = zeros(1, height(points));
    for i = 1:size(meanXYZR, 1)
        xyz = meanXYZR(:, 1:3) - meanXYZR(i, 1:3);
        dist = sum((xyz).^2, 2);
        mask = dist < thresMerge^2;
        mask = mask & maskSize;
        maskPtr = false(1, height(points));
        for j = 1:length(mask)
            if mask(j)
                maskPtr(idx == j) = true;
            end
        end
        % new merged cluster
        if nnz(mask)>0
            count = count + 1;
            idx2(maskPtr) = count;
        end
    end
    idx = idx2;
    
    
    %% score cluster
    clusterSize = zeros(count, 1);
    varXY = zeros(count, 1);
    varR = zeros(count, 1);
    fval = zeros(count, 1);
    meanXYZR = zeros(count, 4);
    for i = 1:count
        mask = idx == i;
        doubleArray = table2array(points(mask, :));
        clusterSize(i) = nnz(mask);
        varXY(i) = sum(var(doubleArray(:, 1:2)));
        varR(i) = var(doubleArray(:, 4));
        fval(i) = mean(abs(doubleArray(:, 5)));
        meanXYZR(i, :) = mean(doubleArray(:, 1:4));
    end
    % temporary fix to small cluster size
    mask = clusterSize > thresClusterSize;
    clusterSize = clusterSize(mask);
    varXY = varXY(mask);
    varR = varR(mask);
    fval = fval(mask);
    meanXYZR = meanXYZR(mask, :);
    fprintf("After clustering:     #clusters = %d\n", length(clusterSize));
    
    score = table;
    score.clusterSize = clusterSize;
    score.varXY = varXY;
    score.varR = varR;
    score.fval = fval;
    score.score = clusterSize./max(clusterSize) - varXY./max(varXY) - varR./max(varR) + fval./max(fval);
    score.meanX = meanXYZR(:, 1);
    score.meanY = meanXYZR(:, 2);
    score.meanZ = meanXYZR(:, 3);
    score.meanR = meanXYZR(:, 4);
    
    %% visualize
    % visualize clusters
    if visPlot == 1

        figure
        colormapExt = [bone(64); parula(64)];
        C = (max(score.meanZ)+1-score.meanZ).*3; % depth: z
        %C = score.clusterSize;
        %C = score.varXY;
        %C = score.varR;
        %C = score.fval;
        %C = score.score;
        
        % outliers
        C90 = quantile(C, 0.9);
        C10 = quantile(C, 0.1);
        C(C < C10) = C10;
        C(C > C90) = C90;
        
        image = image(:, :, 22)+image(:, :, 27)+image(:, :, 33);
        C1 = double(image - min(image, [], 'all'));
        C1 = C1 ./ max(C1) .* 63 + 1;
        h(1) = imshow(C1, [min(C1, [], 'all'), max(C1, [], 'all')]);
        hold on
        C2 = C - min(C);
        C2 = C2 ./ max(C2) .* 63 + 1 + 64;
        h(2) = scatter3(score.meanX, score.meanY, (max(score.meanZ)+1-score.meanZ).*5, 35, C2, 'filled', 'MarkerEdgeColor', [0.8500 0.3250 0.0980]);
        colormap(colormapExt)
        caxis([1, 128]);
        colorbar
    end
    % visualize dots on subplots
    if visPlot == 2

        % number of row, col
        % depth = 12:31; % for full_z=4
        depth = 7:28; % for full_z=3
        ncol = ceil(sqrt(16/9*length(depth)));
        nrow = ceil(length(depth) / ncol);
        % enhance
        t = image(:, :, depth);
        t = t(:);
        if false
            lowerBound = quantile(t, 0.003);
            upperBound = quantile(t, 0.997);
        else
            lowerBound = min(t);
            upperBound = max(t);
        end
        
        % subplot
        figure
        ha = tight_subplot(nrow, ncol);
        count_ = 1;
        for i = depth
            
            axes(ha(count_)); %#ok
            count_ = count_ + 1;
            imshow(image(:, :, i), [lowerBound, upperBound]);
            colormap(pink);
            titleStr = sprintf('z=%d', i);
            title(titleStr);
            hold on
            
            mask = points.z == i;
            points_ = points(mask, :);
            viscircles([points_.x, points_.y], points_.r, 'LineWidth', 0.3, 'Color', [0.2, 0.2, 0.8]);
            % scatter(points_.x, points_.y, 22, 'filled');
        end
    end
    
    %% nested function: filter all eight fields with given mask
    function [] = Filter(mask)
        x = x(mask);
        y = y(mask);
        z = z(mask);
        r = r(mask);
        fval = fval(mask);
        flag = flag(mask);
        startingX = startingX(mask);
        startingY = startingY(mask);
    end
end
