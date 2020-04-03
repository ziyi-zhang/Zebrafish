function [] = SurfXY(testMat)
% by fixing z,r,h, plot the surf of energy w.r.t. x and y

    % specify testMat
    if nargin<1
        testMat = GenerateTestCylinder(1);
    end
    testMat = double(testMat);
    % f_testMat = @(arr)fun(testMat, arr);
    f_testMat = @(arr)GaussianSigmaFun(testMat, arr);
    
    %% correct cylinder attri
    % cyl.x=16; cyl.y=16; cyl.z=5; cyl.r=8; cyl.h=10; % for test image
    % cyl.x=16.5; cyl.y=16.5; cyl.z=17; cyl.r=4.3; cyl.h=7; % for realDot
    cyl.x=45; cyl.y=31; cyl.z=21; cyl.r=4; cyl.h=4; % for 9-marker pt
     
    %%
    delta = 5;
    [xx, yy] = meshgrid(cyl.x-delta:0.3:cyl.x+delta, cyl.y-delta:0.3:cyl.y+delta);
    zz = zeros(size(xx));
    for i = 1:size(zz, 1)
        for j = 1:size(zz, 2)
            zz(j, i) = f_testMat([xx(i, j), yy(i, j), cyl.z, cyl.r, cyl.h]);
        end
    end
    
    % plot
    figure
    mask = zz > 1e5;
    zz(mask) = max(zz(~mask), [], 'all');
    surf(xx, yy, zz, 'FaceAlpha', 0.5);
    colorbar
    title('Energy vs. x and y');
    xlabel('X')
    ylabel('Y')
    zlabel('Energy')
end
