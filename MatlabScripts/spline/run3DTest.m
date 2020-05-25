f = @(x, y, z) x.^2 + x.*z - y.*x;

[x, y, z] = meshgrid(1:20, 1:20, 1:20);
feval = f(x, y, z);

for i = 1:10

    query = rand(1, 3) .* 12 + 3;
    res1 = interp3D(feval, query);
    res2 = interp3(feval, query(2), query(1), query(3), 'spline');
    roundquery = round(query);
    avg = feval(roundquery(1)-1:roundquery(1)+1, roundquery(2)-1:roundquery(2)+1, roundquery(3)-1:roundquery(3)+1);
    avg = mean(avg(:));
    fprintf("avg=%10.3f | b-spline=%10.3f | Interp3=%10.3f\n", avg, res1, res2);
end
