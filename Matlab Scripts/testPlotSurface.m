% plot membrane surface
function [] = testPlotSurface(pt3)
    pt3 = pt3(:, :, 1:30);
    maxNum = max(pt3, [], 3);

    figure
    hold on
    target = zeros(0, 3);
    targetMesh = zeros(size(pt3, 1), size(pt3, 2));
    for i = 1:size(pt3, 1)
        for j = 1:size(pt3, 2)

            numArr = pt3(i, j, :);
            t = quantile(numArr, 0.9);
            for k = 1:size(pt3, 3)

                if (pt3(i, j, k) == maxNum(i, j))
                    target(end+1, :) = [i, j, k];
                    targetMesh(i, j) = k;
                    break;
                end
            end
        end
    end

    %mesh(targetMesh)
    scatter3(target(:, 2), target(:, 1), max(target(:, 3))+1-target(:, 3), 5, max(target(:, 3))+1-target(:, 3), 'filled');
    colormap(parula);
    grid on
    title('Membrane surface')
end
