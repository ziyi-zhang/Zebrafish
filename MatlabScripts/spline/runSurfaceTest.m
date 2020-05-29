f = @(x, y) (x-6).^4 - (y-6).^4 + 2.*x - y;
% f = @(x, y)  x + y;

[x, y] = meshgrid(1:12, 1:12);
z = f(x, y);
z = z + rand(12, 12) .* 3;

% interpolate
[xx, yy] = meshgrid(3:0.2:10, 3:0.2:10);
res = zeros(length(3:0.2:10));
for ix = 1:size(xx, 1)
    for iy = 1:size(xx, 2)

        %{
        t = floor(xx(i));
        xcoef = [basicfunc(xx(i)-t-1), basicfunc(xx(i)-t), basicfunc(xx(i)-t+1), basicfunc(xx(i)-t+2)];
        t = floor(yy(i));
        ycoef = [basicfunc(yy(i)-t-1), basicfunc(yy(i)-t), basicfunc(yy(i)-t+1), basicfunc(yy(i)-t+2)];

        [xcoef_, ycoef_] = meshgrid(xcoef, ycoef);
        %}
        xidx = round((xx(ix, iy)-2.8)/0.2);
        yidx = round((yy(ix, iy)-2.8)/0.2);
        tx = floor(xx(ix, iy));
        ty = floor(yy(ix, iy));
        for j = 1:-1:-2
            for k = 1:-1:-2

                xcoef = basisfunc(xx(ix, iy)-tx+j);
                ycoef = basisfunc(yy(ix, iy)-ty+k);
                res(xidx, yidx) = res(xidx, yidx) + xcoef * ycoef * z(tx-j, ty-k);
            end
        end
    end
end

surf(xx, yy, res,'FaceAlpha',0.5)
hold on
scatter3(x(:), y(:), z(:), 'filled')
