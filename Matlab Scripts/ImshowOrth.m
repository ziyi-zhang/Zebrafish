function [] = ImshowOrth(mat)

    mat = squeeze(mat);
    if size(mat, 1)>size(mat, 2), mat=mat';end
    
    figure
    ha = tight_subplot(3, 1, .01, .01);
    % 1
    axes(ha(1));
    imshow(mat, [min(mat, [], 'all'), max(mat, [], 'all')]);
    colormap pink
    % 2
    axes(ha(2));
    %level = graythresh(mat);
    level = 0.28;
    bw = imbinarize(mat, level);
    imshow(bw);
    % 3
    axes(ha(3));
    bw_sum = sum(bw);
    plot(max(bw_sum)-bw_sum);
    axis tight
end
