% xy accuracy

N = 7;
energy = zeros(N, 1);
xyzrResHist = zeros(N, 4);
xyStartHist = [160, 118; 202, 109; 246, 130; 156, 160; 196, 149; 240, 181; 192, 192];

options = optimoptions(@fminunc,'Display','none','Algorithm','quasi-newton',...
                            'MaxIterations',10,'MaxFunctionEvaluations',800);

img_used = img3;
zRange = 1:16;
N = length(zRange);
for i = 1:7
    energyZHist = zeros(N, 1);
    xResHist = zeros(N, 3);
    for z = zRange
        fsigma = @(arr)GaussianFun(img_used, arr, z);
        [xRes, fval, exitflag, ~] = fminunc(fsigma, [xyStartHist(i, 1), xyStartHist(i, 2), 4], options);
        energyZHist(z) = fval;
        xResHist(z, :) = xRes;
    end
        [minE, idx] = min(energyZHist);
    
        energy(i) = minE;
        xRes = xResHist(idx, :);
        xyzrResHist(i, 1:2) = xRes(1:2);
        xyzrResHist(i, 3) = idx;
        xyzrResHist(i, 4) = xRes(3);
end


