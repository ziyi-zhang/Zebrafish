function [res] = EvaluateCylinder(mat, sampleCyl, samplePeri, weight)
% [res] = EvaluateCylinder(mat, sampleCyl, samplePeri, weight, ratio)
% Calculate the cylindrical evaluation function value with given sample
% points
% 'mat' is a 3d matrix representing a TIFF image
% 'sampleCyl' and 'samplePeri' are the sample points from 'SampleCylinder'
% 'weight' is the weight used by Gaussian Quadrature

    if isempty(sampleCyl)
        % perhaps the cylinder parameters are invalid
        res = double(1e30);
        return;
    end

    evalCyl = interp3(mat, sampleCyl(2, :), sampleCyl(1, :), sampleCyl(3, :), 'spline');
    evalPeri = interp3(mat, samplePeri(2, :), samplePeri(1, :), samplePeri(3, :), 'spline');

    if ~isempty(weight) % 'weight' is only used by Gaussian
        inner = sum(evalCyl .* weight) ./ 2;  % weight here is designed for extended circle, so divide by 2
        extended = sum(evalPeri .* weight);
        res = 2 * inner - extended;
        %res = inner - outer / 0.44;
    else
        % used by 'equadistant'
        res = sum(evalCyl) / length(sampleCyl) - sum(evalPeri) / length(samplePeri);
    end
end
