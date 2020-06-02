% This file is used to experiment the choice of #controls in least square
% b-spline


%controlPts = textread("/Users/ziyizhang/Desktop/NYU/Graphics/Zebrafish/release/debug_50_50_41.log");
controlPts = textread("/Users/ziyizhang/Desktop/NYU/Graphics/Zebrafish/release/debug_30_30_41.log");
feval = pt3;

% params
N_x = size(feval, 1);
N_y = size(feval, 2);
N_z = size(feval, 3);
N = N_x * N_y * N_z;
num_x = 30;
num_y = 30;
num_z = 41;
num = num_x * num_y * num_z;
gap_x = (N_x-1) / (num_x-1);
gap_y = (N_y-1) / (num_y-1);
gap_z = (N_z-1) / (num_z-1);
centers_x = (0:num_x-1).*gap_x + 1; 
centers_y = (0:num_y-1).*gap_y + 1;
centers_z = (0:num_z-1).*gap_z + 1;

% queries
if true
    % xy
    z = 25; % 5 ~ 35
    queries = repmat([0, 0, z], 61*61, 1);
    t = repmat(20:80, 61, 1);
    queries(:, 1) = t(:);
    queries(:, 2) = repmat((20:80)', 61, 1);
    px = 61;
    py = 61;
    
    res2 = feval(20:80, 20:80, z);
    titlestr = "XY with z = " + string(z);
else
    % xz
    y = 30; % 20 ~ 80
    px = 31;
    py = 90;
    queries = repmat([0, y, 0], px*py, 1);
    t = repmat(6:36, py, 1);
    queries(:, 3) = t(:);
    queries(:, 1) = repmat((6:95)', px, 1);
    
    res2 = feval(6:95, y, 6:36);
    res2 = squeeze(res2);
    res2 = res2';
    titlestr = "XZ with y = " + string(y);
end



% visualize
res1 = zeros(px, py);
for i = 1:length(queries)
    
    query = queries(i, :);
    xx = floor((i-1) / py) + 1;
    yy = i - py * (xx - 1);
    
    res1(xx, yy) = interp3DLeastSquare(centers_x, centers_y, centers_z, gap_x, gap_y, gap_z, controlPts, query);
end

figure
set(gcf, 'Position', [100, 100, 1800, 600])
subplot(121)
imagesc(res1);
colorbar
subplot(122)
imagesc(res2);
colorbar

sgtitle(titlestr)
