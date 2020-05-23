function [res, score_] = ScorePoints(inStruct, visPlot, pt1)
% After running 'GlobalImage', score the resultant points

    if nargin<2
        visPlot = false;
    end

    x = inStruct.x;
    y = inStruct.y;
    z = inStruct.z;
    r = inStruct.r;
    fval = inStruct.fval;
    flag = inStruct.flag;
    startingX = inStruct.startingX;
    startingY = inStruct.startingY;
    
    % want fval < -0.05
    mask = fval < -0.1;
    Filter(mask);
    
    % want radius <= 7
    mask = r <= 7;
    Filter(mask);
    
    % Kmeans
    %K = 20;
    K = floor(length(x)/10);
    options = statset('UseParallel', 0);
    [idx, c] = kmeans([x, y], K, 'Options', options, 'replicates', 10);
    
    % merge clusters
    mergedCluster = 1:K;
    newK = K;
    for i = 1:K
        for j = i+1:K
            if (mergedCluster(i) == mergedCluster(j)), continue;end
            if (norm(c(i, :) - c(j, :)) < 1)
                mergedCluster(j) = i;
                newK = newK - 1;
            end
        end
    end
    
    % Score function
    score = zeros(K, 4);
    for i = 1:K
        mask = idx == i;
        xx = x(mask);
        yy = y(mask);
        zz = z(mask);
        rr = r(mask);
        fval_ = fval(mask);
        flag_ = flag(mask);

        clusterSize = length(xx);
        varXY = sum(var([xx, yy]));
        varR = var(rr);
        meanFval = mean(fval_);
        score(i, :) = [clusterSize, varXY, varR, meanFval];
    end
    
    % want varXY < 1
    mask = (score(:, 2) < 1) & (score(:, 1) > 4);
    
    % create a struct
    res = table;
    res.x = x;
    res.y = y;
    res.z = z;
    res.r = r;
    res.fval = fval;
    res.flag = flag;
    res.startingX = startingX;
    res.startingY = startingY;
    % res.cluster = idx(mask);
    % res.clusterX = c(mask, 1);
    % res.clusterY = c(mask, 2);
    %}
    score = score(mask, :);
    c = c(mask, :);
    score_.clusterSize = score(:, 1);
    score_.varXY = score(:, 2);
    score_.varR = score(:, 3);
    score_.meanFval = score(:, 4);

    % final score
    score = score - min(score);
    score = score ./ max(score);
    score_.score = 2 * score(:, 1) - score(:, 2) - score(:, 3) + score(:, 4);
    
    % plot
    if visPlot
        figure
        colormapExt = [bone(64); parula(64)];
        %C = score_.score;
        C = abs(score(:, 4));

        image = pt1(:, :, 22)+pt1(:, :, 27)+pt1(:, :, 33);
        C1 = double(image - min(image, [], 'all'));
        C1 = C1 ./ max(C1) .* 63 + 1;
        h(1) = imshow(C1, [min(C1, [], 'all'), max(C1, [], 'all')]);
        hold on
        C2 = C - min(C);
        C2 = C2 ./ max(C2) .* 63 + 1 + 64;
        h(2) = scatter(c(:, 1), c(:, 2), 35, C2, 'filled', 'MarkerEdgeColor', [0.8500 0.3250 0.0980]);
        colormap(colormapExt)
        caxis([1, 128]);
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
