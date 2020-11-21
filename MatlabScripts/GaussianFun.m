function [res] = GaussianFun(mat, cylArray, z)

    cyl.x = cylArray(1);
    cyl.y = cylArray(2);
    cyl.z = z;  % no longer optimize 'z'
    cyl.r = cylArray(3);
    cyl.h = 2;  % fix height as 2
    
      % Do this once before calling this function:
      % mat = double(mat);
      % mat = mat ./ max(mat, [], 'all');  % normalize to 0-1
    [sampleCyl, samplePeri, weight] = SampleCylinder(cyl, mat, false, 'gaussian');
    
    % evaluate cylinder
    if isempty(sampleCyl)  % ETH data, limit r_max
        % perhaps the cylinder parameters are invalid
        res = double(1);
        return;
    end
    % Create small matrix for interpolation
    padLength = 3;
    xMin = max(1           , floor(min(samplePeri(1, :)))-padLength);
    xMax = min(size(mat, 2), ceil(max(samplePeri(1, :)))+padLength);
    yMin = max(1           , floor(min(samplePeri(2, :)))-padLength);
    yMax = min(size(mat, 1), ceil(max(samplePeri(2, :)))+padLength);
    mat = mat(yMin:yMax, xMin:xMax, :);
    sampleCyl(1, :) = sampleCyl(1, :) - xMin + 1;
    sampleCyl(2, :) = sampleCyl(2, :) - yMin + 1;
    samplePeri(1, :) = samplePeri(1, :) - xMin + 1;
    samplePeri(2, :) = samplePeri(2, :) - yMin + 1;

    evalResCyl = interp3(mat, sampleCyl(1, :), sampleCyl(2, :), sampleCyl(3, :), 'spline');
    evalResPeri = interp3(mat, samplePeri(1, :), samplePeri(2, :), samplePeri(3, :), 'spline');
    
    % res
    alpha = 3;
    K = 3;  % K = R/r, radius of peri over radius of cyl
    % Note: The weight here is QaussianWeight*constants*SubtractionSigma
    res = sum(evalResCyl .* weight) - alpha * sum(evalResPeri .* weight);  % K^2 offsets
    res = -1 * res;  % if the markers are white
end
