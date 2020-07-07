% // interp visualization

figure
hold on
Z = reshape(t_line(:, 3), 90, 24);
[X, Y] = meshgrid(1:90, 1:24);
surf(X-1, 90-Y, -Z')

% set axis xy equal
h = get(gca, 'DataAspectRatio');
set(gca, 'DataAspectRatio', [1 1 1/20])
