% helper function of CPP optimization test
% Interp result vs. [x, y]

fid = fopen('../../release/debug_interpVis_test.log');
% fid = fopen('../../release/debug_interpVis_image.log');
status = 0;
cacheVec = [];


while ~feof(fid)

    line = fgetl(fid);
    
    if ~isempty(line) && line(1)=='>'
        status = 1;
        continue;
    end

    if status == 1
        % read another line
        rawVec = sscanf(line, "%f");
        
        cacheVec = [cacheVec; rawVec'];
    end
end

xmin = min(cacheVec(:, 1));
xmax = max(cacheVec(:, 1));
ymin = min(cacheVec(:, 2));
ymax = max(cacheVec(:, 2));
gap = 1;
[Xq, Yq] = meshgrid(xmin:gap:xmax, ymin:gap:ymax);
[X, Y] = meshgrid(xmin:xmax, ymin:ymax);
Z = reshape(cacheVec(:, 3), ymax-ymin+1, xmax-xmin+1);
Vq = interp2(X, Y, Z, Xq, Yq, 'spline');
figure
hold on
surf(Xq, Yq, -Vq);
grid on

% axis equal
set(gca,'DataAspectRatio',[1 1 1/25])
