Vin = t;
N = size(Vin, 1);
M = N*(N-1)/2;

distRec = zeros(M, 1);
xyRec = zeros(M, 2);
count = 0;
for i = 1:N-1
    for j = i+1:N
        
        count = count + 1;
        distRec(count) = norm(Vin(i, :) - Vin(j, :));
        xyRec(count, :) = [i, j];
    end
end
