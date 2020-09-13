function [] = SubplotDepthCorrection(mat)

    filter = mat == 1;
    mat_without_one = mat(~filter);
    mat(filter) = max(max(mat_without_one));
    [~, minIdx] = min(mat');
    
    % average
    newMat = zeros(size(mat));
    for row = 1:size(mat, 1)
        for col=1:size(mat, 2)
            newMat(row, col) = FindAvg(mat, row, col);
        end
    end
    
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
    end
    sgtitle("Energy for markers optimized with different depths", 'FontSize', 15);
end


function [res] = FindAvg(mat, i, j)

    cols = size(mat, 2);
    t1 = max(1, j-2);
    t2 = min(cols, j+2);
    
    subvec = mat(i, t1:t2);
    res = (sum(subvec) - min(subvec) - max(subvec)) / (t2-t1-1);
end
