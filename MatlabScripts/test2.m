
N = 27;
energy = zeros(N, 1);
xResHist = zeros(N, 3);

options = optimoptions(@fminunc,'Display','none','Algorithm','quasi-newton',...
                            'MaxIterations',10,'MaxFunctionEvaluations',800);
for z = 1:N
    
    fsigma = @(arr)GaussianFun(img, arr, z);
    [xRes, fval, exitflag, ~] = fminunc(fsigma, [162, 118, 4], options);
    
    energy(z) = fval;
    xResHist(z, :) = xRes;
end
