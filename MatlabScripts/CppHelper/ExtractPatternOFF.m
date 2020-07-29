function [] = ExtractPatternOFF(imagePath)
% [] = ExtractPatternOFF(imagePath)
% Read the pattern TIFF image and generate an "OFF" file with extracted pattern 
% location. This is used by the ICP stage in GUI.

    % read the image
    t = Tiff(imagePath, 'r');
    image = t.read();

    % standard circle detection
    BW = imbinarize(image);
    s = regionprops(BW);

       % mask = [s.Area] <= 11; % dont want lines
    centers = [s.Centroid];
    center = [centers(1:2:length(centers)); centers(2:2:length(centers))];
       % center = center(:, mask);
        %mask = (center(1, :) < 400) & (center(2, :) < 400);
        %center = center(:, mask);
    
    if true
        figure
        imshow(image)
        hold on
        scatter(center(1, :), center(2, :));
        title('extracted pattern');
    end
    
    N = size(center, 2);
    loc = zeros(N, 3);
    loc(:, 1:2) = center';

    % write to OFF
    [path, name, ~] = fileparts(imagePath);
    outputName = strcat(name, ".off");
    outputFullName = fullfile(path, outputName);
    fileID = fopen(outputFullName, 'w');
    fprintf(fileID, 'OFF\n');
    fprintf(fileID, '%d %d %d\n\n', N, 0, 0);
    fprintf(fileID, '%d %d %d\n', loc');
    fclose(fileID);
end
