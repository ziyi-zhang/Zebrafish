% convergence zone

N = 49;
energy = zeros(N, 1);
xResHist = zeros(N, 3);
xyStartHist = zeros(N, 2);

options = optimoptions(@fminunc,'Display','none','Algorithm','quasi-newton',...
                            'MaxIterations',10,'MaxFunctionEvaluations',800);

x0 = 160;
y0 = 118;
i = 0;
z = 9;
for x = x0-3:x0+3
    for y = y0-3:y0+3
    
        i = i + 1;
        fsigma = @(arr)GaussianFun(img, arr, z);
        [xRes, fval, exitflag, ~] = fminunc(fsigma, [x, y, 4], options);

        energy(i) = fval;
        xResHist(i, :) = xRes;
        xyStartHist(i, :) = [x, y];
    end
end



figure
hold on
scatter(xResHist(:, 1), xResHist(:, 2))
scatter(xyStartHist(:, 1), xyStartHist(:, 2))
for i = 1:N
    line([xyStartHist(i, 1), xResHist(i, 1)], [xyStartHist(i, 2), xResHist(i, 2)], 'Marker', 'none')
end
