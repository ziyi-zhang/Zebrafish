function [] = SubplotDepthCorrection(mat)

    filter = mat == 1;
    mat_without_one = mat(~filter);
    mat(filter) = max(max(mat_without_one));
    
    % average
    newMat = zeros(size(mat));
    for row = 1:size(mat, 1)
        for col=1:size(mat, 2)
            newMat(row, col) = FindAvg(mat, row, col);
        end
    end
    
    % again
    newMat2 = zeros(size(mat));
    for row = 1:size(mat, 1)
        for col=1:size(mat, 2)
            newMat2(row, col) = FindAvg(newMat, row, col);
        end
    end
    

    % critical line
    [~, minIdx] = min(newMat2');
    minIdxLeft = minIdx - 8;
    minIdxRight = minIdx + 8;
    minIdxLeft(minIdxLeft < 1) = 1;
    minIdxRight(minIdxRight > size(mat, 2)) = size(mat, 2);
    
    figure
    N = size(mat, 1);
    rows = floor(sqrt(N)*3/4);
    cols = ceil(N / rows);
    [ha, ~] = tight_subplot(rows, cols);
    for i = 1:N
        axes(ha(i));
        plot(1:size(mat, 2), mat(i, :));
        hold on
        plot(1:size(mat, 2), newMat(i, :))
        plot(1:size(mat, 2), newMat2(i, :))
        % two vertical lines
        plot([minIdxLeft(i), minIdxLeft(i)], [min(newMat(i, :)), max(newMat(i, :))], '--', 'Color', [0.9 0.1 0.1])
        plot([minIdxRight(i), minIdxRight(i)], [min(newMat(i, :)), max(newMat(i, :))], '--', 'Color', [0.9 0.1 0.1])
    end
    sgtitle("Energy for markers optimized with different depths", 'FontSize', 15);
    
    figure
    N = size(mat, 1);
    rows = floor(sqrt(N)*3/4);
    cols = ceil(N / rows);
    [ha, ~] = tight_subplot(rows, cols);
    for i = 1:N
        axes(ha(i));
        
        ttt = conv(newMat2(i, :), [-1, 2, -1], 'valid');
        plot(2:size(newMat2, 2)-1, ttt);
        ylimVec = [-6e-3, 6e-3];
        ylim(ylimVec);
        hold on
        % two vertical lines
        plot([minIdxLeft(i), minIdxLeft(i)], ylimVec, '--', 'Color', [0.9 0.1 0.1])
        plot([minIdxRight(i), minIdxRight(i)], ylimVec, '--', 'Color', [0.9 0.1 0.1])
    end
    sgtitle("2nd Energy for markers optimized with different depths", 'FontSize', 15);
    %
end


function [res] = FindAvg(mat, i, j)

    cols = size(mat, 2);
    t1 = max(1, j-2);
    t2 = min(cols, j+2);
    
    subvec = mat(i, t1:t2);
    res = (sum(subvec) - min(subvec) - max(subvec)) / (t2-t1-1);
end
