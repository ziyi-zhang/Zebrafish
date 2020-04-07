
% kmeans test
searchRange = 10:40;
varArray = zeros(size(searchRange));
count = 1;
for i = searchRange

    [idx, c, sumd] = kmeans([xx yy], i);
    varArray(count) = sum(sumd);
    count = count + 1;
end
