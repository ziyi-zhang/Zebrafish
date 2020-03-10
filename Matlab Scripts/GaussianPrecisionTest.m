function [] = GaussianPrecisionTest()
% This function is used to compare the precision of Gaussian Quadrature and
% equadistant integral

    f = @(x, y) sqrt(2-x.^2-y.^2);
    % f = @(x, y) x.^2 + y.^2;

    %% Gaussian - with proposal 1 to calc subtraction
    [xyCyl, xyPeri, weightCyl, weightPeri] = SampleCylinderHelper('gaussian');
    ratio = 1.3;
    GaussianResult = ratio^2 * (sum(f(xyCyl(1, :), xyCyl(2, :)) .* weightCyl) - sum(f(xyPeri(1, :), xyPeri(2, :)) .* weightPeri));
    %GaussianResult = ratio^2 * (sum(f(xyCyl(1, :), xyCyl(2, :)) .* weightCyl) + sum(f(xyPeri(1, :), xyPeri(2, :)) .* weightPeri));
    GaussianEvals = length(xyCyl)+length(xyPeri);
    
    %% Gaussian - with proposal 2 to calc subtraction
    OuterResult = ratio^2 * (sum(f(xyCyl(1, :), xyCyl(2, :)) .* weightCyl) + sum(f(xyPeri(1, :), xyPeri(2, :)) .* weightPeri));
    InnerResult = sum(f(xyCyl(1, :)./ratio, xyCyl(2, :)./ratio) .* weightCyl) + sum(f(xyPeri(1, :)./ratio, xyPeri(2, :)./ratio) .* weightPeri);
    Gaussian2Result = 2 * InnerResult - OuterResult;
    fprintf("Gaussian-2 result   %.10f with %d evaluations.\n", Gaussian2Result, 2*GaussianEvals);
    
    %% Equadistant
    [xyCyl, xyPeri] = SampleCylinderHelper('equadistant', 200, 60);
    EquaResult = sum(f(xyCyl(1, :), xyCyl(2, :))) - sum(f(xyPeri(1, :), xyPeri(2, :)));
    %EquaResult = sum(f(xyCyl(1, :), xyCyl(2, :))) + sum(f(xyPeri(1, :), xyPeri(2, :)));
    EquaEvals = length(xyCyl)+length(xyPeri);
    EquaResult = EquaResult ./ EquaEvals * ratio^2*pi;
    
    %% Log
    
    fprintf("Gaussian-1 result:  %.10f with %d evaluations.\n", GaussianResult, GaussianEvals);
    fprintf("Equadistant result: %.10f with %d evaluations.\n", EquaResult, EquaEvals);
end
