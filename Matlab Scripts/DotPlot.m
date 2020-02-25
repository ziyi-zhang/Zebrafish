function [] = DotPlot(mat)
    
    min_ = min(mat, [], 'all');
    max_ = max(mat, [], 'all');
    pink_ = pink(double(max_-min_+1));
    
    figure
    hold on
    for i = 1:4:size(mat, 1)
        for j = 1:4:size(mat, 2)
            for k = 1:size(mat, 3)
            
                scatter3(i, j, k, 4, pink_(mat(i, j, k)-min_+1));
            end
        end
    end
end
