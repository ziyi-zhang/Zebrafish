function [p] = GetRawLoc(cdata)

    BW = imbinarize(cdata(:, :, 1));
    s = regionprops(BW);

    mask = [s.Area] <= 11; % dont want lines
    centers = [s.Centroid];
    center = [centers(1:2:length(centers)); centers(2:2:length(centers))];
    center = center(:, mask);
    
    if true
        figure
        imshow(cdata)
        hold on
        scatter(center(1, :), center(2, :));
        title('Get dot location from raw pattern');
    end
    
    p = zeros(3, size(center, 2));
    p(1:2, :) = center;
end
