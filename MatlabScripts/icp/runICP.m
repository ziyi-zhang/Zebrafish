% function [] = runICP(p)

q = GetRawLoc(cdata);
p = [score.meanX'; score.meanY'; zeros(1, height(score))];
% p = p .* 1.84;
ratio = 1.7:0.002:1.9;

minErr = 1e10;
bestScale = 0;
errorArray = zeros(size(ratio));
for i = 1:length(ratio)

    pp = p .* ratio(i);
    [TR, TT, ER] = icp(q, pp);
    errorArray(i) = ER(end);
    if (ER(end) < minErr)
        minErr = ER(end);
        bestScale = ratio(i);
    end
end
fprintf('Cost: %f\n', minErr);
fprintf('Scale: %f\n', bestScale);
p = p .* bestScale;
[TR, TT, ER] = icp(q, p);

visResult = true;
if visResult

    figure
    hold on
    % scatter3(p(1, :), p(2, :), p(3, :));
    scatter3(q(1, :), q(2, :), q(3, :), 'MarkerEdgeColor', [0 0.4470 0.7410]);
    pp = TR * p + TT;
    scatter3(pp(1, :), pp(2, :), pp(3, :), 'MarkerEdgeColor', [0.8500 0.3250 0.0980]);
    legend('Raw Pattern', 'Calculated Results');
end
