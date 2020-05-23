function [] = TestMatScript(testMat)
% Check correctness on test matrix

    % specify testMat
    if nargin<1
        testMat = GenerateTestCylinder(1);
    end
    testMat = double(testMat);
    f_testMat = @(arr)fun(testMat, arr);
    fsigma_testMat = @(arr)GaussianSigmaFun(testMat, arr);
    %% correct cylinder attri
    % cyl.x=16; cyl.y=16; cyl.z=5; cyl.r=8; cyl.h=10;  % for test image
    % cyl.x=17; cyl.y=17; cyl.z=17; cyl.r=4; cyl.h=7;
    cyl.x=45; cyl.y=31; cyl.z=21; cyl.r=4; cyl.h=3; % for 9-marker pt
    
    %% f_testMat
    tic
    %
    f1 = figure;
    % x
    f1(1) = subplot(2, 2, 1);
    range = cyl.x-4:0.2:cyl.x+4;
    res = zeros(size(range));
    count = 1;
    for i = range
        res(count) = f_testMat([i, cyl.y, cyl.z, cyl.r, cyl.h]);
        count = count+1;
    end
    plot(range, res);
    hold on
    yl = ylim;
    plot([cyl.x, cyl.x], yl, '-.', 'Color', '#D95319');
    grid on
    title('Energy vs cyl.x');
    % z
    f1(2) = subplot(2, 2, 2);
    range = cyl.z-4:0.2:cyl.z+4;
    res = zeros(size(range));
    count = 1;
    for i = range
        res(count) = f_testMat([cyl.x, cyl.y, i, cyl.r, cyl.h]);
        count = count+1;
    end
    plot(range, res);
    hold on
    yl = ylim;
    plot([cyl.z, cyl.z], yl, '-.', 'Color', '#D95319');
    grid on
    title('Energy vs cyl.z');
    % r
    f1(3) = subplot(2, 2, 3);
    range = 1:0.2:cyl.r+3;
    res = zeros(size(range));
    count = 1;
    for i = range
        res(count) = f_testMat([cyl.x, cyl.y, cyl.z, i, cyl.h]);
        count = count+1;
    end
    plot(range, res);
    hold on
    yl = ylim;
    plot([cyl.r, cyl.r], yl, '-.', 'Color', '#D95319');
    grid on
    title('Energy vs cyl.r');
    % h
    f1(4) = subplot(2, 2, 4);
    t = max(1, cyl.h-4);
    range = t:0.2:cyl.h+4;
    res = zeros(size(range));
    count = 1;
    for i = range
        res(count) = f_testMat([cyl.x, cyl.y, cyl.z, cyl.r, i]);
        count = count+1;
    end
    plot(range, res);
    hold on
    yl = ylim;
    plot([cyl.h, cyl.h], yl, '-.', 'Color', '#D95319');
    grid on
    title('Energy vs cyl.h');
    sgtitle('Subtraction by two evaluations - f-testMat');
    %
    toc
    
    %% fsigma_testMat
    tic
    f2 = figure;
    % x
    f2(1) = subplot(2, 2, 1);
    range = cyl.x-4:0.2:cyl.x+4;
    res = zeros(size(range));
    count = 1;
    for i = range
        res(count) = fsigma_testMat([i, cyl.y, cyl.z, cyl.r, cyl.h]);
        count = count+1;
    end
    plot(range, res);
    hold on
    yl = ylim;
    plot([cyl.x, cyl.x], yl, '-.', 'Color', '#D95319');
    grid on
    title('Energy vs cyl.x');
    % z
    f2(2) = subplot(2, 2, 2);
    range = cyl.z-4:0.2:cyl.z+4;
    res = zeros(size(range));
    count = 1;
    for i = range
        res(count) = fsigma_testMat([cyl.x, cyl.y, i, cyl.r, cyl.h]);
        count = count+1;
    end
    plot(range, res);
    hold on
    yl = ylim;
    plot([cyl.z, cyl.z], yl, '-.', 'Color', '#D95319');
    grid on
    title('Energy vs cyl.z');
    % r
    f2(3) = subplot(2, 2, 3);
    range = 1:0.2:cyl.r+3;
    res = zeros(size(range));
    count = 1;
    for i = range
        res(count) = fsigma_testMat([cyl.x, cyl.y, cyl.z, i, cyl.h]);
        count = count+1;
    end
    plot(range, res);
    hold on
    yl = ylim;
    plot([cyl.r, cyl.r], yl, '-.', 'Color', '#D95319');
    grid on
    title('Energy vs cyl.r');
    % h
    f2(4) = subplot(2, 2, 4);
    t = max(1, cyl.h-4);
    range = t:0.2:cyl.h+4;
    res = zeros(size(range));
    count = 1;
    for i = range
        res(count) = fsigma_testMat([cyl.x, cyl.y, cyl.z, cyl.r, i]);
        count = count+1;
    end
    plot(range, res);
    hold on
    yl = ylim;
    plot([cyl.h, cyl.h], yl, '-.', 'Color', '#D95319');
    grid on
    title('Energy vs cyl.h');
    sgtitle('Subtraction function with one evaluaion - fsigma-testMat');
    toc
end
