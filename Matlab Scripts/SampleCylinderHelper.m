function [xyCyl, xyPeri] = SampleCylinderHelper(n, e, visualization)
% [xyCyl, xyPeri] = SampleCylinderHelper(n, e, visualization)
% Helper function to generate equadistant sample points for 'SampleCylinder'
% 'n' interior layers, excluding the origin
% 'e' exterior layers
% 'visualization' is a flag default as false

    if nargin<3, visualization=false;end

    xyCyl = [0, 0];
    for i = 1:n
    
        theta = 1/i;
        m = ceil(2*pi/theta);
        theta = 2*pi/m;
        for j = 0:m-1
            xyCyl = [xyCyl; i*cos(theta*j), i*sin(theta*j)]; %#ok
        end
    end
    xyCyl = xyCyl .* (2) ./ (2*n+1);
    xyCyl = transpose(xyCyl);
    xyPeri = [];
    for i = n+1:n+e
    
        theta = 1/i;
        m = ceil(2*pi/theta);
        theta = 2*pi/m;
        for j = 0:m-1
            xyPeri = [xyPeri; i*cos(theta*j), i*sin(theta*j)]; %#ok
        end
    end
    xyPeri = xyPeri .* (2) ./ (2*n+1);
    xyPeri = transpose(xyPeri);

    if visualization
        figure
        scatter(xyCyl(1, :), xyCyl(2, :));
        hold on
        scatter(xyPeri(1, :), xyPeri(2, :));
        axis equal
    end
end
