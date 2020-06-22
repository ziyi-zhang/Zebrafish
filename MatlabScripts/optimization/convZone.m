% helper function of CPP optimization test
% convergence zone

fid = fopen('../../release/debug_deg2_convZone.log');
status = 0;
cacheVec = [];

figure
hold on
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
        cacheVec = cacheVec([1, size(cacheVec, 1)], :);
        plot(cacheVec(:, 2), cacheVec(:, 3), 'LineStyle', '-', 'Marker', '.');
        cacheVec = [];
        status = 0;
    end
end
