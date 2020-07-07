% helper function of CPP optimization test
% plot convergence zone

function [] = convZone()
% fid = fopen('../../release/debug_deg2_quad3_test.log');
% fid = fopen('../../release/debug_deg2_quad3_pt1.log');

% fid = fopen('../../build/debug.log');
fid = fopen('../../release/debug.log');

status = 0;
count = 0;
correctCnt = 0;
correctCenter = [12.36, 11.38];
cacheVec = [];
evalVec = [];

figure
hold on
grid on
while ~feof(fid)

    line = fgetl(fid);
    
    if status == 1
        rawVec = sscanf(line, "%f");
        cacheVec = [cacheVec; rawVec'];
    end
    
    % status
    if ~isempty(line) && line(1)=='>'
        status = 1;
    end
    if ~isempty(line) && line(1)=='<'
        status = -1;
    end

    % plot
    if status == -1
    
        if isempty(cacheVec), continue;end
        
        count = count + 1;
        if (norm(correctCenter - cacheVec(end, 2:3)) < 0.1), correctCnt = correctCnt + 1;end
        evalVec(count) = size(cacheVec, 1);

        cacheVec = cacheVec([1, size(cacheVec, 1)], :);
        plot(cacheVec(:, 2), cacheVec(:, 3), 'LineStyle', '-', 'Marker', '.');
        cacheVec = [];
        status = 0;
    end
end

fprintf("Starting point count = %d\n", count);
fprintf("Converged count = %d\n", correctCnt);
end