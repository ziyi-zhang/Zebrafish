function [res] = EvaluateCylinder(mat, sampleCyl, samplePeri, weight)
% [res] = EvaluateCylinder(mat, sampleCyl, samplePeri, weight)
% Calculate the cylindrical evaluation function value with given sample
% points
% 'mat' is a 3d matrix representing a TIFF image
% 'sampleCyl' and 'samplePeri' are the sample points from 'SampleCylinder'
% 'weight' is the weight used by Gaussian Quadrature

    if isempty(sampleCyl)
        % perhaps the cylinder parameters are invalid
        res = realmax;
        return;
    end

    evalCyl = interp3(mat, sampleCyl(2, :), sampleCyl(1, :), sampleCyl(3, :), 'spline');
    evalPeri = interp3(mat, samplePeri(2, :), samplePeri(1, :), samplePeri(3, :), 'spline');

    if ~isempty(weight)
        inner = sum(evalCyl .* weight);
        outer = sum(evalPeri .* weight) - inner;
        res = inner - outer / 0.44;
    else
        res = sum(evalCyl) - sum(evalPeri) * (length(sampleCyl) / length(samplePeri));
    end
end
