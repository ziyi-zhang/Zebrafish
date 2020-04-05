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
        res = double(10);
        return;
    end
    evalRes = interp3(mat, sample(1, :), sample(2, :), sample(3, :), 'spline');
    
    % res
    res = sum(evalRes .* weight); % Note: The weight here is QaussianWeight*constants*SubtractionSigma
end
