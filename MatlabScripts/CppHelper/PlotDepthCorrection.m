function [] = PlotDepthCorrection(mat)

    filter = mat == 1;
    mat_without_one = mat(~filter);
    mat(filter) = max(max(mat_without_one));
    [~, minIdx] = min(mat');
    
    figure
    imagesc(mat);
    colorbar
    hold on
    scatter(minIdx, 1:size(mat, 1), 'red', 'o', 'filled');
    title("Energy for markers optimized with different depths", 'FontSize', 15);
end
