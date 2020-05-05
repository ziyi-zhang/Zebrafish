function [res] = runOpticalFlow(image1, image2, points)
% Run HS3D with previous frame 'image1' and subsequet frame 'image2'
% 'points' is a N-by-3 matrix representing the markers in 'image1'
% 'res' is a N-by-3 matrix representing new locations in 'image2'

    addpath './OpticalFlow/';
    smoothness = 1;
    iterations = 50;
    
    [ux, uy, uz] = HS3D(image1, image2, smoothness, iterations);
    rounded = round(points);
    % boundary
    rounded(rounded == 0) = 1;
    rounded(:, 3) = rounded(:, 3) + 2;  % to middle layer of the cylinder
    mask = rounded(:, 3) > size(image1, 3);
    rounded(mask, 3) = size(image1, 3);
    % res
    res = points;
    for i = 1:size(points, 1)
    
       % 'points' should be xyz
       %{
       correction = [0, 0, 0];
       correction(1) = interp3(ux, rounded(i, 2), rounded(i, 1), rounded(i, 3), 'linear');
       correction(2) = interp3(uy, rounded(i, 2), rounded(i, 1), rounded(i, 3), 'linear');
       correction(3) = interp3(uz, rounded(i, 2), rounded(i, 1), rounded(i, 3), 'linear');
       %}
       %{
       for dx = -1:1:1
           for dy = -1:1:1
               correction = correction + [ux(rounded(i, 2)+dx, rounded(i, 1)+dy, rounded(i, 3)),...
                                          uy(rounded(i, 2)+dx, rounded(i, 1)+dy, rounded(i, 3)),...
                                          uz(rounded(i, 2)+dx, rounded(i, 1)+dy, rounded(i, 3))]; % center
           end
       end
       correction = correction ./ 9.0;
       %}
       correction = [ux(rounded(i, 2), rounded(i, 1), rounded(i, 3)),...
                      uy(rounded(i, 2), rounded(i, 1), rounded(i, 3)),...
                      uz(rounded(i, 2), rounded(i, 1), rounded(i, 3))];
       res(i, :) = res(i, :) + correction;
    end
end
