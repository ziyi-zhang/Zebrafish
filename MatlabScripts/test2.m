
minIdxArray = zeros(size(p, 2), 1);
distArray = zeros(size(p, 2), 1);

for i = 1:size(p, 2)
    
    minDistSq = 1e8;
    minIdx = 0;
    for j = 1:size(q, 2)
    
        % distSq = norm( p(:, i) - q(:, j) );
        distSq = norm( [p(1, i)-q(2, j), p(2, i)-q(1, j), p(3, i)-q(3, j)] );
        if (distSq < minDistSq)
            minDistSq = distSq;
            minIdx = j;
        end
    end
    
    minIdxArray(i) = minIdx;
    distArray(i) = minDistSq;
end
