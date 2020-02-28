function [resCyl, resPeri] = SampleCylinder(cyl, n)
% [resCyl, resPeri] = SampleCylinder(cyl, n)
% 'cyl' is a struct representing a cylinder
% 'n' denotes the level of sample
% 'resCyl' is of size [3 by n^3] representing n^3 points in cylinder
% 'resPeri' is of size [3 by n^3] representing n^3 peripheral points
% Sample cylinder:
% cyl.x = 10, cyl.y = 20, cyl.z = 30, cyl.r = 5, cyl.h = 40;

    if nargin<2, n=10;end

    x = cyl.x;
    y = cyl.y;
    z = cyl.z;
    r = cyl.r;
    h = cyl.h;

    % depth array
    pad = h/n/2;
    zArray = linspace(z-pad, z+h+pad, n+2);
    zArray = zArray(2:end-1);

    % xy (rc) array
    xyCylArray = zeros(2, n*n);
    xyPeriArray = zeros(2, n*n);
    %{
    count = 0;
    for i = 0:n-1 % theta

        theta = 2*pi/n * i;
        delta = [sin(theta); cos(theta)] .* r ./ n; % note x is row, not col
        for j = 1:n % rho

            count = count + 1;
            xyCylArray(:, count) = [x; y] + (j-0.5) .* delta;
            xyPeriArray(:, count) = [x; y] + (j+n-0.5) .* delta;
        end
    end
    %}
    % a faster version of the above commented block
    for i = 0:n-1 % theta
        
        theta = 2*pi/n * i;
        delta = [sin(theta); cos(theta)] .* r ./ n;
        xyCylArray(:, i*n+1:(i+1)*n) = repmat([x; y], 1, n) + ((1:n)-0.5) .* delta;
        xyPeriArray(:, i*n+1:(i+1)*n) = repmat([x; y], 1, n) + (((1+n):(2*n))-0.5) .* delta;
    end

    % res
    zArray = repmat(zArray, n*n, 1);
    zArray = transpose(zArray(:));
    resCyl = repmat(xyCylArray, 1, n);
    resCyl = [resCyl; zArray];
    resPeri = repmat(xyPeriArray, 1, n);
    resPeri = [resPeri; zArray];

    % debug
    debug = false;
    if debug
        figure
        hold on
        scatter3(resCyl(2, :), resCyl(1, :), resCyl(3, :));
        scatter3(resPeri(2, :), resPeri(1, :), resPeri(3, :));
        axis tight
    end
end
