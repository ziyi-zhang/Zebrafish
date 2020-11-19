function [cluster, points] = FilterPoints(inMat, image, visPlot)
% After running cpp version code, score the resultant points
    
    if nargin<3
        visPlot = 0;
    end

    % filter config
    fvalThres = -0.05;
    radiusThres = 7;
    iterMax = 15;
    membraneThres = 0.7; % binary threshold
    HistPeriThres = 0.85; % quantile threshold (of peripheral samples' variance)
    % cluster config
    thresDist = 0.01;
    thresClusterSize = 25;
    
    % retrive data
    x = inMat(:, 1);
    y = inMat(:, 2);
    z = inMat(:, 3);
    r = inMat(:, 4);
    fval = inMat(:, 5);
    iter = inMat(:, 6);
    startingX = inMat(:, 7);
    startingY = inMat(:, 8);
    fprintf("Raw #starting points            = %d\n", length(x));
    
    % want fval < fvalThres
    mask = fval < fvalThres;
    Filter(mask);
    fprintf("After filter fval:      #points = %d\n", length(x));
    
    % want radius <= radiusThres
    mask = r <= radiusThres;
    Filter(mask);
    fprintf("After filter radius:    #points = %d\n", length(x));
    
    % want iter <= iterMax
    mask = iter < iterMax;
    Filter(mask);
    fprintf("After iterL #points = %d\n", length(x));
    
    
    % want peri color bright (in membrane area)
    %{
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
    %}
    
    % hist of peripheral area
    %{
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
    %}
    

    % create a struct
    points = table;
    points.x = x;
    points.y = y;
    points.z = z;
    points.r = r;
    points.fval = fval;
    points.iter = iter;
    points.startingX = startingX;
    points.startingY = startingY;
    
    %% cluster
    [~, sortIdx] = sort(x);
    clusterSet = DJSet(length(x));
    for i = 1:length(x)
        for j = i+1:length(x)
            
            if (x(sortIdx(j)) - x(sortIdx(i)) > thresDist), break;end

            dist_square = (x(sortIdx(j)) - x(sortIdx(i)))^2 + (y(sortIdx(j)) - y(sortIdx(i)))^2;
            if (dist_square < thresDist^2)
                clusterSet.union(sortIdx(i), sortIdx(j));
            end
        end
    end
    
    belongIdx = zeros(length(x), 1);
    for i = 1:length(x)
        belongIdx(i) = clusterSet.find(i);
    end
    clusterIdx = unique(belongIdx);
    numClusters = length(clusterIdx);
    fprintf("#clusters = %d\n", numClusters);
    clusterMat = zeros(numClusters, 7);
    numClustersSize = 0;
    for i = 1:numClusters
        t = clusterIdx(i);
        mask = belongIdx == t;
        if (nnz(mask) > thresClusterSize)
            numClustersSize = numClustersSize + 1;
        else
            continue;
        end
        
        clusterMat(numClustersSize, 1) = mean(x(mask));
        clusterMat(numClustersSize, 2) = mean(y(mask));
        clusterMat(numClustersSize, 3) = mean(z(mask));
        clusterMat(numClustersSize, 4) = mean(r(mask));
        clusterMat(numClustersSize, 5) = nnz(mask);
        clusterMat(numClustersSize, 6) = mean(fval(mask));
        clusterMat(numClustersSize, 7) = mean(iter(mask));
    end
    fprintf("After size filter #clusters = %d\n", numClustersSize);

    %% Prepare to output
    cluster = table;
    cluster.meanX = clusterMat(1:numClustersSize, 1);
    cluster.meanY = clusterMat(1:numClustersSize, 2);
    cluster.meanZ = clusterMat(1:numClustersSize, 3);
    cluster.meanR = clusterMat(1:numClustersSize, 4);
    cluster.count = clusterMat(1:numClustersSize, 5);
    cluster.energy = clusterMat(1:numClustersSize, 6);
    cluster.meanIter = clusterMat(1:numClustersSize, 7);
    
    %% Visualization
    % visualize clusters
    if visPlot == 1

        figure
        colormapExt = [bone(64); parula(64)];
        C = (max(cluster.meanZ)+1-cluster.meanZ).*3; % depth: z
        
        % outliers
        C90 = quantile(C, 0.9);
        C10 = quantile(C, 0.1);
        C(C < C10) = C10;
        C(C > C90) = C90;
        
        image = image(:, :, 23)+image(:, :, 30)+image(:, :, 34)+image(:, :, 38);
        % image = image(:, :, 22)+image(:, :, 27)+image(:, :, 33);
        C1 = double(image - min(image, [], 'all'));
        C1 = C1 ./ max(C1) .* 63 + 1;
        h(1) = imshow(C1, [min(C1, [], 'all'), max(C1, [], 'all')]);
        hold on
        C2 = C - min(C);
        C2 = C2 ./ max(C2) .* 63 + 1 + 64;
        h(2) = scatter3(cluster.meanY+1, cluster.meanX+1, (max(cluster.meanZ)+1-cluster.meanZ).*5, ...
                        35, C2, 'filled', 'MarkerEdgeColor', [0.8500 0.3250 0.0980]);
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
        iter = iter(mask);
        startingX = startingX(mask);
        startingY = startingY(mask);
    end
end
