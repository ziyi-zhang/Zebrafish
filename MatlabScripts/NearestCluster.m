function [score] = NearestCluster(points, image)

    %% cluster
    thresDist = 0.15;
    thresClusterSize = 10;
    
    count = 0;
    idx = zeros(1, height(points));
    for i = 1:height(points)
    
        if (idx(i) > 0), continue;end
        xy = table2array(points(:, 1:2));
        xy = xy - xy(i, :);
        dist = xy(:, 1).^2 + xy(:, 2).^2;
        mask = dist < thresDist;
        if (nnz(mask) < thresClusterSize), continue;end
        % new cluster
        count = count + 1;
        idx(mask) = count;
    end
    
    %% score
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
    if false
    
        figure
        colormapExt = [bone(64); parula(64)];
        C = score.clusterSize;
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
        h(2) = scatter(score.meanX, score.meanY, 35, C2, 'filled', 'MarkerEdgeColor', [0.8500 0.3250 0.0980]);
        colormap(colormapExt)
        caxis([1, 128]);
    end
    % visualize dots on subplots
    if true
        
        % number of row, col
        depth = 12:31;
        ncol = ceil(sqrt(16/9*length(depth)));
        nrow = ceil(length(depth) / ncol);
        % enhance
        t = image(:, :, depth);
        t = t(:);
        lowerBound = quantile(t, 0.003);
        upperBound = quantile(t, 0.997);
        
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
            viscircles([points_.x, points_.y], points_.r, 'LineWidth', 0.5, 'Color', [0.2, 0.2, 0.8]);
            % scatter(points_.x, points_.y, 22, 'filled');
        end
    end
end
