function [res] = GaussianSigmaFun(mat, cylArray)

    cyl.x = cylArray(1);
    cyl.y = cylArray(2);
    cyl.z = cylArray(3);
    cyl.r = cylArray(4);
    cyl.h = cylArray(5);
    mat = mat ./ max(mat, [], 'all');  % normalize to 0-1
    [sample, ~, weight] = SampleCylinder(cyl, mat, false, 'gaussiansigma');
    
    % evaluate cylinder
    if isempty(sample)
        % perhaps the cylinder parameters are invalid
        res = double(1e30);
        return;
    end
    evalRes = interp3(mat, sample(2, :), sample(1, :), sample(3, :), 'spline');
    
    % res
    res = sum(evalRes .* weight); % Note: The weight here is QaussianWeight*constants*SubtractionSigma
end
