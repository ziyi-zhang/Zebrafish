% helper function of CPP optimization test
% energy function vs. [x, y]

fid = fopen('../../release/debug.log');
cacheVec = [];


while ~feof(fid)

    line = fgetl(fid);
    
    if ~isempty(line) && line(1)=='>'
        
        % read another line
        line = fgetl(fid);
        rawVec = sscanf(line, "%f");
        
        cacheVec = [cacheVec; rawVec(2), rawVec(3), rawVec(5)];
    end
end

xmin = min(cacheVec(:, 1));
xmax = max(cacheVec(:, 1));
ymin = min(cacheVec(:, 2));
ymax = max(cacheVec(:, 2));
[Xq, Yq] = meshgrid(xmin:0.3:xmax, ymin:0.3:ymax);
[X, Y] = meshgrid(12-4:12+4, 11-4:11+4);
Z = reshape(cacheVec(:, 3), 9, 9);
Vq = interp2(X, Y, Z, Xq, Yq, 'spline');
figure
hold on
surf(Xq, Yq, Vq);
