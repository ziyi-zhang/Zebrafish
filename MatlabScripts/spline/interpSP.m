function [res] = interpSP(y, query)
% Test cubic b-spline interpolation

    res = zeros(size(query));
    for i = 1:length(query)
    
        t = query(i);
        k = floor(t);
        % k-1
        tt = t - k + 1;
        res(i) = res(i) + (- 1/6 * tt^3 + tt^2 -2*tt + 4/3) * y(k-1);
        % k
        tt = t - k;
        res(i) = res(i) + (1/2*tt^3 - tt^2 + 2/3) * y(k);
        % k+1
        tt = t - k - 1;
        res(i) = res(i) + (-1/2*tt^3 - tt^2 + 2/3) * y(k+1);
        % k+2
        tt = t - k - 2;
        res(i) = res(i) + (1/6*tt^3 + tt^2 + 2*tt + 4/3) * y(k+2);
    end
end

% match the result of spmak
