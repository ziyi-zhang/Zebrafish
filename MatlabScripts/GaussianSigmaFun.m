function [res] = GaussianSigmaFun(mat, cylArray, z)

    cyl.x = cylArray(1);
    cyl.y = cylArray(2);
    cyl.z = z;  % no longer optimize 'z'
    cyl.r = cylArray(3);
    cyl.h = 4;  % fix height as 4
    
      % Do this once before calling this function:
      % mat = double(mat);
      % mat = mat ./ max(mat, [], 'all');  % normalize to 0-1
    [sample, ~, weight] = SampleCylinder(cyl, mat, false, 'gaussiansigma');
    
    % evaluate cylinder
    if isempty(sample)
        % perhaps the cylinder parameters are invalid
        res = double(1);
        return;
    end
    % Create small matrix for interpolation
    padLength = 2;
    xMin = max(1           , floor(min(sample(1, :)))-padLength);
    xMax = min(size(mat, 2), ceil(max(sample(1, :)))+padLength);
    yMin = max(1           , floor(min(sample(2, :)))-padLength);
    yMax = min(size(mat, 1), ceil(max(sample(2, :)))+padLength);
    mat = mat(yMin:yMax, xMin:xMax, :);
    sample(1, :) = sample(1, :) - xMin + 1;
    sample(2, :) = sample(2, :) - yMin + 1;
    evalRes = interp3(mat, sample(1, :), sample(2, :), sample(3, :), 'spline');
    
    % res
    res = sum(evalRes .* weight);  % Note: The weight here is QaussianWeight*constants*SubtractionSigma
    res = -1 * res;  % if the markers are white
end
