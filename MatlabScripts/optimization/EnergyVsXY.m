% helper function of CPP optimization test
% Energy function vs. [x, y]

% fid = fopen('../../release/debug_cylinderVis_test.log');
% fid = fopen('../../release/debug_cylinderVis_pt1.log');

% fid = fopen('../../build/debug.log');
fid = fopen('../../release/debug.log');

status = 0;
cacheVec = [];

figure
hold on
grid on

subPltCnt = 1;
radiusArr = [2.0, 2.8, 3.6, 4.4];
while ~feof(fid)

    line = fgetl(fid);
    
    % status
    if ~isempty(line) && line(1)=='>'
        status = 1;
        continue;
    end
    if ~isempty(line) && line(1)=='<'
        status = -1;
    end
       
    if status == 1
        % read another line
        rawVec = sscanf(line, "%f");
        
        cacheVec = [cacheVec; rawVec'];
    end
    
    if status == -1
    
        % plot
        subplot(1, 4, subPltCnt);
        xmin = min(cacheVec(:, 1));
        xmax = max(cacheVec(:, 1));
        ymin = min(cacheVec(:, 2));
        ymax = max(cacheVec(:, 2));
        gap = 0.8;
        [Xq, Yq] = meshgrid(xmin:gap:xmax, ymin:gap:ymax);
        [X, Y] = meshgrid(xmin:xmax, ymin:ymax);
        Z = reshape(cacheVec(:, 3), ymax-ymin+1, xmax-xmin+1);
        Vq = interp2(X, Y, Z, Xq, Yq, 'spline');
        surf(Xq, Yq, Vq, 'FaceAlpha', 0.9);
        view(2);
        % axis equal
        set(gca,'DataAspectRatio',[1 1 1/30])
        axis([3, 33, 3, 76]);
        %axis([3, 28, 3, 28]);
        titleStr = sprintf("Radius = %f", radiusArr(subPltCnt));
        title(titleStr);
        
        status = 0;
        subPltCnt = subPltCnt + 1;
        cacheVec = [];
    end
end

sgtitle("Energy vs. XY with fixed radius");
