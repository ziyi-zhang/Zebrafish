% detect all targets
% [centers, radii] = step2(image);

function [centers, radii] = step2(image, radiusLow, radiusHigh)

    if nargin <= 1
        radiusLow = 2;
        radiusHigh = 6;
    end
    frames = length(image);
    centers = cell(frames, 1);
    radii = cell(frames, 1);

    for f = 1:frames
        [center, radius] = imfindcircles(max(image{f}, [], 3), [radiusLow, radiusHigh], 'Sensitivity', 0.95);
        
        gap = 15;
        N = size(center, 1);
        mask = true(N, 1);
        for i = 1:N
            
            if (center(i, 1) <= gap) || (center(i, 2) <= gap) || (center(i, 1)+gap>size(image{f}, 1)) || (center(i, 2)+gap>size(image{f}, 2))
                mask(i) = false;
            end
        end
        centers{f} = center(mask, :);
        radii{f} = radius(mask, :);
    end
end
