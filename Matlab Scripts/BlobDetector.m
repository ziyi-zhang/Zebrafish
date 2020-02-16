function [centers, radii] = BlobDetector(image)

    BW_threshold = 0.99;
    radii_min = 2.5;
    radii_max = 6.5;

    image_bw = image > quantile(image(:), BW_threshold);
    % regionprop
    stats = regionprops('table', image_bw, 'Centroid', 'MajorAxisLength', 'MinorAxisLength');
    centers = stats.Centroid;
    radii = mean([stats.MajorAxisLength stats.MinorAxisLength], 2) ./ 2;
    % mask
    mask = (radii > radii_min) & (radii < radii_max);
    centers = centers(mask, :);
    radii = radii(mask);
end
