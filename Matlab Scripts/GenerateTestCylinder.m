function [testMat] = GenerateTestCylinder(visualize)
% testMat = GenerateTestCylinder(1);
% Generate a 3D matrix with several cylinders 
% DEBUG USE

    if nargin<1, visualize=false;end
    matSize = [32, 32, 41];
    testMat = ones(matSize(1), matSize(2), matSize(3));
    
    % 1 perfect cylinder
    cyl.x=16; cyl.y=16; cyl.z=5; cyl.r=8; cyl.h=10;
    for i = cyl.x-cyl.r:cyl.x+cyl.r
        for j = cyl.y-cyl.r:cyl.y+cyl.r
            if norm([i, j]-[cyl.x, cyl.y])<=cyl.r
                testMat(i, j, cyl.z:cyl.z+cyl.h) = 0;
            end
        end
    end
    
    % 2 truncated cone
    %{
    cyl.x=200; cyl.y=500; cyl.z=5; cyl.r=20; cyl.h=10;
    for i = cyl.x-cyl.r:cyl.x+cyl.r
        for j = cyl.y-cyl.r:cyl.y+cyl.r
            for k = cyl.z:cyl.z+cyl.h
                if norm([i, j]-[cyl.x, cyl.y])<cyl.r-(k-cyl.z)
                    testMat(i, j, k) = 0;
                end
            end
        end
    end
    %}
    
    
    % visualize
    if visualize
        figure
        title('Visualization of test cylinders');
        hold on
        target = [];
        for i = 1:matSize(1)
            for j = 1:matSize(2)
                for k = 1:matSize(3)
                    if testMat(i, j, k)==0
                        target = [target; i, j, k]; %#ok
                    end
                end
            end
        end
        scatter3(target(:, 2), target(:, 1), target(:, 3), 8, 'filled');
        hold on
        viscircles([cyl.y, cyl.x], cyl.r, 'LineWidth', 0.5)
        xlim([1 matSize(1)])
        ylim([1 matSize(2)])
        zlim([-1 matSize(3)])
        view(3)
        grid on
    end
end
