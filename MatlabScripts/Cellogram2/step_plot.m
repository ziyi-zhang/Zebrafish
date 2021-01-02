% plot

function [] = step_plot(xyzreHist)
    
    % need a figure
    hold on
    grid on
    frames = length(xyzreHist);
    
    for f = 1:frames
        scatter3(xyzreHist{f}(:, 1), xyzreHist{f}(:, 2), xyzreHist{f}(:, 3));
    end
end
