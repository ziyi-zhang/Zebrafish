function [xHist, inputHist, fvalHist, fvalOldHist] = StartPointEffect(fun)
% Used to plot a figure to see how stable the optimization method is

    options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton',...
                            'MaxIterations',1000,'MaxFunctionEvaluations',1000,...
                            'OptimalityTolerance',1e-10);
    
    figure
    hold on
    N = 100;
    xHist = zeros(N, 5);
    inputHist = zeros(100, 2);
    fvalHist = zeros(N, 5);
    fvalOldHist = zeros(N, 5);
    count = 0;
    for i = 8:2:26
        for j = 8:2:26
            [x, fval, exitflag, output] = fminunc(fun, [i, j, 15, 4, 7], options);
            fval_old = fun([i, j, 15, 4, 7]);
            plot3([j x(2)], [i x(1)], [fval_old fval], '-o');
            
            count = count + 1;
            xHist(count, :) = x;
            inputHist(count, :) = [i, j];
            fvalHist(count, :) = fval;
            fvalOldHist(count, :) = fval_old;
            fprintf('%d %d\n', i, j);
        end
    end
end
