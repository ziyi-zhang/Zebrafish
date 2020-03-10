function [] = GaussianPrecisionTest()

    % f = @(x, y) sqrt(2-x.^2-y.^2);
    f = @(x, y) x.^2 + y.^2;

    %% Gaussian
    [xyCyl, xyPeri, weightCyl, weightPeri] = SampleCylinderHelper('gaussian');
    ratio = 1.3;
    %GaussianResult = sum(f(xyCyl(1, :), xyCyl(2, :)) .* weightCyl) - sum(f(xyPeri(1, :), xyPeri(2, :)) .* weightPeri);
    GaussianResult = sum(f(xyCyl(1, :), xyCyl(2, :)) .* weightCyl) + sum(f(xyPeri(1, :), xyPeri(2, :)) .* weightPeri);
    GaussianEvals = length(xyCyl)+length(xyPeri);
    
    %% Equadistant
    [xyCyl, xyPeri] = SampleCylinderHelper('equadistant', 200, 60);
    %EquaResult = sum(f(xyCyl(1, :), xyCyl(2, :))) - sum(f(xyPeri(1, :), xyPeri(2, :)));
    EquaResult = sum(f(xyCyl(1, :), xyCyl(2, :))) + sum(f(xyPeri(1, :), xyPeri(2, :)));
    EquaEvals = length(xyCyl)+length(xyPeri);
    EquaResult = EquaResult ./ EquaEvals;
    
    %% Log
    
    fprintf("Gaussian result:    %10.7f with %d evaluations.\n", GaussianResult, GaussianEvals);
    fprintf("Equadistant result: %10.7f with %d evaluations.\n", EquaResult, EquaEvals);
end
