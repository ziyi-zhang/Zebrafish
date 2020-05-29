function [res, centers, y] = interpSP(y, query, num)
% Test cubic b-spline interpolation

    if size(y, 1)==1, y=y';end
    
    % least square fitting
    if num > 0
        N = length(y);
        A = zeros(N, num);
        gap = (N-1) / (num-1);
        centers = (0:num-1).*gap + 1;
        for j = 1:num
            center = centers(j);
            for i = 1:N
                
                t = (i-center)/gap;
                A(i, j) = basisfunc(t);
            end
        end
        
        y = linsolve(A'*A, A'*y);
    end
    % Evaluate
    res = zeros(size(query));
    for i = 1:length(query)
    
        t = query(i);
        idx = floor((t-1)/gap) + 1;

            res(i) = res(i) + basisfunc((t-centers(idx-1))/gap) * y(idx-1); 
            res(i) = res(i) + basisfunc((t-centers(idx  ))/gap) * y(idx  ); 
            res(i) = res(i) + basisfunc((t-centers(idx+1))/gap) * y(idx+1); 
            res(i) = res(i) + basisfunc((t-centers(idx+2))/gap) * y(idx+2); 
        %{
        for j = 1:num
            center = centers(j);
            res(i) = res(i) + basisfunc((t-center)/gap) * y(j); 
        end
        %}
    end
end

% match the result of spmak
