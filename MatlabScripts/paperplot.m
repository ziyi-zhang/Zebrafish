% plot for paper

function [] = paperplot(mat)

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
    
    figure
    hold on
    ax = gca;
    ax.FontSize = 15; 
    plot(1:size(mat, 2), mat(1, :));
    plot(1:size(mat, 2), newMat2(1, :))
    xlabel('Z Axis', 'FontSize', 17)
    ylabel('Energy', 'FontSize', 17)
end


function [res] = FindAvg(mat, i, j)

    cols = size(mat, 2);
    t1 = max(1, j-2);
    t2 = min(cols, j+2);
    
    subvec = mat(i, t1:t2);
    res = (sum(subvec) - min(subvec) - max(subvec)) / (t2-t1-1);
end
