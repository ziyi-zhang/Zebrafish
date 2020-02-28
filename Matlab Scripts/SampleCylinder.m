function [resCyl, resPeri] = SampleCylinder(cyl, n, mat, visualize)
% [resCyl, resPeri] = SampleCylinder(cyl, n, mat, visualize)
% 'cyl' is a struct representing a cylinder
% 'n' denotes the level of sample
% 'mat' is the 3D matrix of image. It will only be used to do boundary
% check and can be omitted.
% 'visualize' is a debug flag of whether to plot the sample points
% 'resCyl' is of size [3 by n^3] representing n^3 points in cylinder
% 'resPeri' is of size [3 by n^3] representing n^3 peripheral points
% Sample cylinder:
% cyl.x = 10; cyl.y = 20; cyl.z = 30; cyl.r = 5; cyl.h = 40;

    if nargin<4, visualize=false;end
    if nargin<3, mat=[];end
    if nargin<2, n=10;end
    if n<5, warning('n should be larger than 4');return;end

    x = cyl.x;
    y = cyl.y;
    z = cyl.z;
    r = cyl.r;
    h = cyl.h;
    D = 3; % # difference of interior and exterior points in radial direction

    % range check
    if ~isempty(mat)
        if x < 1 || x > size(mat, 1) || ...
           y < 1 || y > size(mat, 2) || ...
           z < 1 || z > size(mat, 3) || ...
           r < 0 || r > size(mat, 1)/2 || ...
           h < 0 || z+h > size(mat, 3)

            resCyl = [];
            resPeri = [];
            return;
        end
    end
    
    % depth array
    pad = h/n/2;
    zNumber = ceil(h) + 2; % number of depth
    zArray = linspace(z-pad, z+h+pad, zNumber+2);
    zArray = zArray(2:end-1);

    % xy (rc) array
    xyCylArray = zeros(2, n*n);
    xyPeriArray = zeros(2, n*(n-D));
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
        xyPeriArray(:, i*(n-D)+1:(i+1)*(n-D)) = repmat([x; y], 1, n-D) + (((1+n):(2*n-D))-0.5) .* delta;
    end

    % res
    zArrayPeri = repmat(zArray, n*(n-D), 1);
    zArrayPeri = transpose(zArrayPeri(:));
    zArray = repmat(zArray, n*n, 1);
    zArray = transpose(zArray(:));
    resCyl = repmat(xyCylArray, 1, zNumber);
    resCyl = [resCyl; zArray];
    resPeri = repmat(xyPeriArray, 1, zNumber);
    resPeri = [resPeri; zArrayPeri];

    % debug
    if visualize
        figure
        hold on
        scatter3(resCyl(2, :), resCyl(1, :), resCyl(3, :));
        scatter3(resPeri(2, :), resPeri(1, :), resPeri(3, :));
        axis tight
    end
end
