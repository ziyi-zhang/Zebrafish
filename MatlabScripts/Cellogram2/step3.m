% optimize for xyz
% [xyzreHist] = step3(image, centers);

function [res] = step3(image, centers, zResolution)

    if nargin <= 2
        zResolution = 0.5;
    end

    options = optimoptions(@fminunc,'Display','none','Algorithm','quasi-newton',...
                            'MaxIterations',10,'MaxFunctionEvaluations',800);
    frames = length(image);
    slices = size(image{1}, 3);
    res = cell(frames, 1);
    
    tic
    parfor f = 1:frames

        N = size(centers{f}, 1);  % different frames may have different #targets
        xyzreHist = zeros(N, 5);
        
        % optimization
        sumOptimalZIdx1 = floor(slices/2);
        sumOptimalZIdx2 = length(1:zResolution:2);
        
        for k = 1:N
            %% round 1 (rough zRange)
            zRange = 1:slices;
            M = length(zRange);
            energyHist = zeros(M, 1);
            xyrHist = zeros(M, 3);
            % round 1 - guess trial
            guess = floor(sumOptimalZIdx1 / (k));
            fsigma = @(arr)GaussianFun(image{f}, arr, zRange(guess));
            [xRes, fval, ~, ~] = fminunc(fsigma, [centers{f}(k, 1), centers{f}(k, 2), 3], options);
            energyHist(guess) = fval;
            xyrHist(guess, :) = xRes;
            % round 1 - left to right
            cntSubsequentIncrease = 0;
            for i = guess+1:M
                z = zRange(i);
                fsigma = @(arr)GaussianFun(image{f}, arr, z);
                [xRes, fval, ~, ~] = fminunc(fsigma, [centers{f}(k, 1), centers{f}(k, 2), 3], options);
                energyHist(i) = fval;
                xyrHist(i, :) = xRes;
                
                % small optimization
                if fval > energyHist(i-1)
                    cntSubsequentIncrease = cntSubsequentIncrease + 1;
                else
                    cntSubsequentIncrease = 0;
                end
                if (cntSubsequentIncrease == 4)
                    energyHist(i+1:end) = 1.1;
                    break;
                end
            end
            % round 1 - right to left
            cntSubsequentDecrease = 0;
            for i = guess-1:-1:1
                z = zRange(i);
                fsigma = @(arr)GaussianFun(image{f}, arr, z);
                [xRes, fval, ~, ~] = fminunc(fsigma, [centers{f}(k, 1), centers{f}(k, 2), 3], options);
                energyHist(i) = fval;
                xyrHist(i, :) = xRes;
                
                % small optimization
                if fval > energyHist(i+1)
                    cntSubsequentDecrease = cntSubsequentDecrease + 1;
                else
                    cntSubsequentDecrease = 0;
                end
                if (cntSubsequentDecrease == 4)
                    energyHist(1:i-1) = 1.1;
                    break;
                end
            end

            %% find min energy in round 1
            [~, minEIdx] = min(energyHist);
            sumOptimalZIdx1 = sumOptimalZIdx1 + minEIdx;
            x = xyrHist(minEIdx, 1);
            y = xyrHist(minEIdx, 2);
            r = xyrHist(minEIdx, 3);
            idx1 = minEIdx - 1;
            idx2 = minEIdx + 1;
            if (minEIdx == 1), idx1 = 1;end
            if (minEIdx == M), idx2 = slices;end

            %% round 2 (fine zRange)
            zRange2 = idx1:zResolution:idx2;
            M = length(zRange2);
            energyHist2 = zeros(M, 1);
            xyrHist2 = zeros(M, 3);
            for i = 1:M
                z = zRange2(i);
                if z == zRange(idx1)
                    energyHist2(i) = energyHist(idx1);
                    xyrHist2(i, :) = xyrHist(idx1, :);
                    continue;
                elseif z == zRange(idx2)
                    energyHist2(i) = energyHist(idx2);
                    xyrHist2(i, :) = xyrHist(idx2, :);
                    continue;
                elseif z == zRange(minEIdx)
                    energyHist2(i) = energyHist(minEIdx);
                    xyrHist2(i, :) = xyrHist(minEIdx, :);
                    continue;
                end
                fsigma = @(arr)GaussianFun(image{f}, arr, z);
                [xRes, fval, ~, ~] = fminunc(fsigma, [x, y, r], options);
                energyHist2(i) = fval;
                xyrHist2(i, :) = xRes;
            end

            %% find min energy in round 2
            [minE, minEIdx] = min(energyHist2);
            sumOptimalZIdx2 = sumOptimalZIdx2 + minEIdx;
            x = xyrHist2(minEIdx, 1);
            y = xyrHist2(minEIdx, 2);
            r = xyrHist2(minEIdx, 3);
            xyzreHist(k, 1) = x;
            xyzreHist(k, 2) = y;
            xyzreHist(k, 3) = zRange2(minEIdx);
            xyzreHist(k, 4) = r;
            xyzreHist(k, 5) = minE;
            
            fprintf("frame=%d target=%d:   E=%.3f R=%.1f\n", f, k, minE, r);
        end

        res{f} = xyzreHist;
    end
    toc
end
