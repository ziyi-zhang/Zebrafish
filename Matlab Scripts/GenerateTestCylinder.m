function [testMat] = GenerateTestCylinder(visualize)
% Generate a 3D matrix with several cylinders 
% DEBUG USE

    if nargin<1, visualize=false;end
    testMat = ones(1024, 1024, 41);
    
    % 1
    cyl.x=200; cyl.y=300; cyl.z=5; cyl.r=20; cyl.h=10;
    for i = cyl.x-cyl.r:cyl.x+cyl.r
        for j = cyl.y-cyl.r:cyl.y+cyl.r
            if norm([i, j]-[cyl.x, cyl.y])<cyl.r
                testMat(i, j, cyl.z:cyl.z+cyl.h) = 0;
            end
        end
    end
    
    % visualize
    if visualize
        figure
        title('Visualization of test cylinders');
        hold on
        target = [];
        for i = 1:1024
            for j = 1:1024
                for k = 1:41
                    if testMat(i, j, k)==0
                        target = [target; i, j, k];
                    end
                end
            end
        end
        scatter3(target(:, 2), target(:, 1), target(:, 3), 8, 'filled');
        xlim([1 1024])
        ylim([1 1024])
        zlim([1 41])
        view(3)
        grid on
    end
end
