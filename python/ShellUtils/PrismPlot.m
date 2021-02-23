function [] = PrismPlot(pr)

    figure
    hold on
    grid on

    connect = [1, 2; 2, 3; 3, 1; 1, 4; 2, 5; 3, 6; 4, 5; 5, 6; 6, 4];
    for i = 1:3
        plot3([pr(connect(i, 1), 1), pr(connect(i, 2), 1)], ...
              [pr(connect(i, 1), 2), pr(connect(i, 2), 2)], ...
              [pr(connect(i, 1), 3), pr(connect(i, 2), 3)], '-', 'Color', 'b');
    end
    for i = 7:9
        plot3([pr(connect(i, 1), 1), pr(connect(i, 2), 1)], ...
              [pr(connect(i, 1), 2), pr(connect(i, 2), 2)], ...
              [pr(connect(i, 1), 3), pr(connect(i, 2), 3)], '-', 'Color', 'r');
    end
    for i = 4:6
        plot3([pr(connect(i, 1), 1), pr(connect(i, 2), 1)], ...
              [pr(connect(i, 1), 2), pr(connect(i, 2), 2)], ...
              [pr(connect(i, 1), 3), pr(connect(i, 2), 3)], '-');
    end
    
    extra = [3, 5];
    for i = 1:1
        plot3([pr(extra(i, 1), 1), pr(extra(i, 2), 1)], ...
              [pr(extra(i, 1), 2), pr(extra(i, 2), 2)], ...
              [pr(extra(i, 1), 3), pr(extra(i, 2), 3)], '-');
    end
end
