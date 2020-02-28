function [res] = EvaluateCylinder(mat, sampleCyl, samplePeri)
% [res] = EvaluateCylinder(mat, sampleCyl, samplePeri)
% Calculate the cylindrical evaluation function value with given sample
% points
% 'mat' is a 3d matrix representing a TIFF image
% 'sampleCyl' and 'samplePeri' are the sample points from 'SampleCylinder'

    if isempty(sampleCyl)
        % perhaps the cylinder parameters are invalid
        res = intmax;
        return;
    end
    
    evalCyl = interp3(mat, sampleCyl(2, :), sampleCyl(1, :), sampleCyl(3, :), 'linear');
    evalPeri = interp3(mat, samplePeri(2, :), samplePeri(1, :), samplePeri(3, :), 'linear');

    res = sum(evalCyl) - sum(evalPeri);
end
