function [res, image] = BlobDetector(image_log, showImage, visBlob)

    BW_threshold = 0.993;
    radii_min = 2.5;
    radii_max = 6.5;

    if nargin<3, visBlob=false;end
    if nargin<2, showImage=false;end
    N = size(image_log, 3);
    res = cell(2, N);
    for i = 1:N
        image = image_log(:, :, i);
        image_bw = image > quantile(image(:), BW_threshold);
        % regionprop
        stats = regionprops('table', image_bw, 'Centroid', 'MajorAxisLength', 'MinorAxisLength');
        centers = stats.Centroid;
        radii = mean([stats.MajorAxisLength stats.MinorAxisLength], 2) ./ 2;
        % mask
        mask = (radii > radii_min) & (radii < radii_max);
        res{1, i} = centers(mask, :);
        res{2, i} = radii(mask);
    end
    
    if showImage
        ImshowTiff(image_log(:, :, 1));
        viscircles(res{1, 1}, res{2, 1});
    end
    
    if visBlob
        image = zeros(size(image_log(:, :, 1)));
        allCenters = [];
        for i = 1:N
            allCenters = cat(1, allCenters, res{1, i});
        end
        allCenters = allCenters';

        for i = 1:size(image, 1)
            for j = 1:size(image, 2)
                image(j, i) = nnz(vecnorm(allCenters - [i; j]) < 3);
            end
        end

        %{
        for i = 1:size(allCenters, 2)
            r = round(allCenters(2, i));
            c = round(allCenters(1, i));
            image(r, c) = image(r, c) + 1;
        end
        %}
        %{
        mask = image < 1;
        image(mask) = 0;
        
        figure
        imshow(image, [min(image, [], 'all'), max(image, [], 'all')]);
        colormap pink
        %}
    end
end
