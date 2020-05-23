function [res] = ClusterResults(inRes, mat)

    xHist = inRes.xHist;
    fvalHist = inRes.fvalHist;
    
    mask = fvalHist < -0.05;
    xHist = xHist(mask, :);
    
    res = [];
    thresDist = 0.8;
    thresCount = 6;
    countMat = zeros(size(mat, 1), size(mat, 2));
    for y = 1:size(countMat, 1)
        for x = 1:size(countMat, 2)
        
            diffVec = xHist(:, 1:2) - repmat([x, y], size(xHist, 1), 1);
            distArr = diffVec(:, 1).^2 + diffVec(:, 2).^2;
            count = nnz(distArr < thresDist^2);
            countMat(y, x) = count;
            if count >= thresCount
                res = [res; x, y];
            end
        end
    end
end
