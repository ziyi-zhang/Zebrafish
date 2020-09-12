function [] = SubplotDepthCorrection(mat)

    filter = mat == 1;
    mat_without_one = mat(~filter);
    mat(filter) = max(max(mat_without_one));
    [~, minIdx] = min(mat');
    
    figure
    N = size(mat, 1);
    rows = floor(sqrt(N)*3/4);
    cols = ceil(N / rows);
    [ha, ~] = tight_subplot(rows, cols);
    for i = 1:N
        axes(ha(i));
        plot(1:size(mat, 2), mat(i, :));
    end
    sgtitle("Energy for markers optimized with different depths", 'FontSize', 15);
end
