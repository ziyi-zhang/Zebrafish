f = @(x, y) (x-10).^3 - (y-10).^3 - 2.*x.^2 + y.^2;
% f = @(x, y)  x + y;

[x, y] = meshgrid(1:20, 1:20);
z = f(x, y);
% z = z + rand(12, 12) .* 3;

% least square fitting
N_x = size(z, 1);
N_y = size(z, 2);
N = N_x * N_y;
num_x = 20;
num_y = 20;
num = num_x * num_y;
gap_x = (N_x-1) / (num_x-1);
gap_y = (N_y-1) / (num_y-1);
centers_x = (0:num_x-1).*gap_x + 1;
centers_y = (0:num_y-1).*gap_y + 1;
A = zeros(N, num);
for jx = 1:num_x
    for jy = 1:num_y
        j = (jx-1)*num_y + jy;
        center_x = centers_x(jx);
        center_y = centers_y(jy);
        for ix = 1:N_x
            for iy = 1:N_y
                i = (ix-1)*N_y + iy;
                
                A(i, j) = basisfunc((ix-center_x)/gap_x) * basisfunc((iy-center_y)/gap_y);
            end
        end
    end
end
inputPts = z'; 
inputPts = inputPts (:);
controlPts = linsolve(A'*A, A'*inputPts);

% interpolate
[xx, yy] = meshgrid(3:0.2:18, 3:0.2:18);
res = zeros(length(3:0.2:18));
for ix = 1:size(xx, 1)
    for iy = 1:size(xx, 2)

        %{
        t = floor(xx(i));
        xcoef = [basicfunc(xx(i)-t-1), basicfunc(xx(i)-t), basicfunc(xx(i)-t+1), basicfunc(xx(i)-t+2)];
        t = floor(yy(i));
        ycoef = [basicfunc(yy(i)-t-1), basicfunc(yy(i)-t), basicfunc(yy(i)-t+1), basicfunc(yy(i)-t+2)];

        [xcoef_, ycoef_] = meshgrid(xcoef, ycoef);
        %}
        xres = round((xx(ix, iy)-2.8)/0.2);
        yres = round((yy(ix, iy)-2.8)/0.2);
        idx_x = floor( (xx(ix, iy)-1)/gap_x ) + 1;
        idx_y = floor( (yy(ix, iy)-1)/gap_y ) + 1;
        for j = -1:2
            for k = -1:2

                my_idx_x = idx_x + j;
                my_idx_y = idx_y + k;
                xcoef = basisfunc( (xx(ix, iy)-centers_x(my_idx_x)) / gap_x );
                ycoef = basisfunc( (yy(ix, iy)-centers_y(my_idx_y)) / gap_y );
                res(xres, yres) = res(xres, yres) + xcoef * ycoef * controlPts( (my_idx_x-1) * num_y + my_idx_y );
            end
        end
    end
end

surf(xx, yy, res,'FaceAlpha',0.4)
hold on
scatter3(x(:), y(:), z(:), 'filled')
[x_control, y_control] = meshgrid(centers_x, centers_y);
x_control = x_control';
y_control = y_control';
scatter3(x_control(:), y_control(:), controlPts, 'filled');
xlim([2.5 18.5])
ylim([2.5 18.5])
zlim([-600 400])
set(gcf, 'Position',  [100, 100, 1200, 1000])
legend("interpolation/fitting surface", "Input points", "Control points")
