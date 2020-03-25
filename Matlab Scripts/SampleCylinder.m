function [resCyl, resPeri, weight] = SampleCylinder(cyl, mat, visualize, method)
% [resCyl, resPeri, weightCyl, weightPeri] = SampleCylinder(cyl, mat, visualize, method)
% 'cyl' is a struct representing a cylinder
% 'mat' is the 3D matrix of image. It will only be used to do boundary
% check and is optional.
% 'visualize' is a debug flag of whether to plot the sample points
% 'resCyl' is of size [3 by n^3] representing n^3 points in cylinder
% 'resPeri' is of size [3 by n^3] representing n^3 peripheral points
% 'weight' is the weight used by Gaussian Quadrature
% Sample cylinder:
% cyl.x = 10; cyl.y = 20; cyl.z = 30; cyl.r = 5; cyl.h = 10;


    persistent xyArray xyArraySigma subtractionArray 
    if nargin<4, method='gaussiansigma';end
    if nargin<3, visualize=false;end
    if nargin<2, mat=[];end

    x = cyl.x;
    y = cyl.y;
    z = cyl.z;
    r = cyl.r;
    h = cyl.h;

    %% Quassian Quadrature with subtraction function sigma
    if strcmp(method, 'gaussiansigma')
        H = 10; % desired height layers
        
        % range check
        if ~isempty(mat)
            if x < r || x > size(mat, 1)-r || ...
               y < r || y > size(mat, 2)-r || ...
               z < 1 || z > size(mat, 3)-1 || ...
               r < 1 || r > 19 || ...
               h < 1 || z+h > size(mat, 3)

                resCyl = [];
                resPeri = [];
                weight = [];
                return;
            end
        end
        
        % depth array
        pad = h / H / 2;
        zArray = linspace(z-pad, z+h+pad, H+2);
        zArray = zArray(2:end-1);
        
        % xy (rc) array
        if isempty(xyArraySigma)
            [xyArraySigma.cyl, ~, xyArraySigma.weight, ~] = SampleCylinderHelper('gaussiansigma');
        end
        if isempty(subtractionArray)
            subtractionArray = SubtractionSigma(xyArraySigma.cyl);
        end
        xyCylArray = xyArraySigma.cyl .* r;
        xyCylArray = xyCylArray + repmat([x; y], 1, size(xyCylArray, 2));
        
        % res
        zCylArray = repmat(zArray, size(xyCylArray, 2), 1);
        zCylArray = transpose(zCylArray(:));
        resCyl = repmat(xyCylArray, 1, length(zArray));
        resCyl = [resCyl; zCylArray];
        weight = repmat(xyArraySigma.weight, 1, length(zArray));
        % weight = weight .* (r^2) .* h;  % Note: The constants are multipied here
        weight = weight .* repmat(subtractionArray, 1, length(zArray)); % subtraction array multiplied here
        weight = weight .* (h ./ H);  % equadistant 1D integral (z-dim)
        weight = weight .* 2 .* r^2;  % gaussian disk jocabian for extended circle (not density jacobian)
        weight = weight ./ (r^2*h*pi);  % Vperi
        resPeri = [];
    end

    %% Gaussian Quadrature (does not use smooth subtraction)
    if strcmp(method, 'gaussian')
        H = 10;
        % H = ceil(h) + 1; % desired height layers
        
        % range check
        if ~isempty(mat)
            if x < r || x > size(mat, 1)-r || ...
               y < r || y > size(mat, 2)-r || ...
               z < 1 || z > size(mat, 3)-1 || ...
               r < 1 || r > 19 || ...
               h < 1 || z+h > size(mat, 3)

                resCyl = [];
                resPeri = [];
                weight = [];
                return;
            end
        end

        % depth array
        pad = h / H / 2;
        zArray = linspace(z-pad, z+h+pad, H+2);
        zArray = zArray(2:end-1);
        
        % xy (rc) array
        if isempty(xyArray)
            [xyArray.cyl, xyArray.peri, xyArray.weight, ~] = SampleCylinderHelper('gaussian');
        end
        xyCylArray = xyArray.cyl .* r;
        xyCylArray = xyCylArray + repmat([x; y], 1, size(xyCylArray, 2));
        xyPeriArray = xyArray.peri .* r;
        xyPeriArray = xyPeriArray + repmat([x; y], 1, size(xyPeriArray, 2));
        
        % res
        zCylArray = repmat(zArray, size(xyCylArray, 2), 1);
        zCylArray = transpose(zCylArray(:));
        resCyl = repmat(xyCylArray, 1, length(zArray));
        resCyl = [resCyl; zCylArray];
        zPeriArray = repmat(zArray, size(xyPeriArray, 2), 1);
        zPeriArray = transpose(zPeriArray(:));
        resPeri = repmat(xyPeriArray, 1, length(zArray));
        resPeri = [resPeri; zPeriArray];
        weight = repmat(xyArray.weight, 1, length(zArray));
        weight = weight * (h ./ H);  % equadistant 1D integral (z-dim)
        weight = weight .* 2 .* (r^2);  % gaussian disk jocabian for extended circle (not density jacobian)
        weight = weight ./ (r^2*h*pi);  % Vperi
        % weight = weight ./ h;
    end

    %% equadistant
    if strcmp(method, 'equadistant')
        if isempty(xyArray), xyArray = cell(50);end
        R = 10; % desired interior layers
        H = ceil(h) + 1; % desired height layers
        E = 4; % layer of exterior sample points

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
        pad = h / H / 2;
        zArray = linspace(z-pad, z+h+pad, H+2);
        zArray = zArray(2:end-1);

        % xy (rc) array
        if isempty(xyArray{R})
            [xyArray{R}.cyl, xyArray{R}.peri] = SampleCylinderHelper('equadistant', R, E);
        end
        xyCylArray = xyArray{R}.cyl .* r;
        xyCylArray = xyCylArray + repmat([x; y], 1, size(xyCylArray, 2));
        xyPeriArray = xyArray{R}.peri .* r;
        xyPeriArray = xyPeriArray + repmat([x; y], 1, size(xyPeriArray, 2));

        % res
        zCylArray = repmat(zArray, size(xyCylArray, 2), 1);
        zCylArray = transpose(zCylArray(:));
        resCyl = repmat(xyCylArray, 1, length(zArray));
        resCyl = [resCyl; zCylArray];
        zPeriArray = repmat(zArray, 1, size(xyPeriArray, 2));
        zPeriArray = transpose(zPeriArray(:));
        resPeri = repmat(xyPeriArray, 1, length(zArray));
        resPeri = [resPeri; zPeriArray];
        weight = [];
    end

    %% polar
    if strcmp(method, 'polar')
        n = 8; % the level of sample
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
        weight = [];
    end

    %% debug
    if visualize
        figure
        hold on
        colorTable = [0.8, 0.8, 0.2; 0.8, 0.2, 0.8];
        filter = weight > 0;
        absWeight = abs(weight);
        %scatter3(resCyl(2, :), resCyl(1, :), resCyl(3, :), ceil(absWeight*200), colorTable(filter+1, :));
        scatter3(resCyl(2, :), resCyl(1, :), resCyl(3, :));
        scatter3(resPeri(2, :), resPeri(1, :), resPeri(3, :));
        axis tight
        view(3)
    end
end
