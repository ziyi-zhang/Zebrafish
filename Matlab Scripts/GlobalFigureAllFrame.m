function [quasiPointsHist, flowPointsHist, imageFrame] = GlobalFigureAllFrame(imageStruct, points, visualize)
% After getting the locations of the first frame, use optical flow and
% quasi-newton to get the movements of cylinders
% 'imageStruct' should be a struct containing metadata with only one
% channel
% 'points' are the clusters returned by 'ScorePoints2'
% Example [quasiPointsHist, flowPointsHist, imageFrame] = GlobalFigureAllFrame(ImageC1, [score.meanX, score.meanY, score.meanZ]);

    if nargin<3
        visualize = false;
    end

    % split frames
    numFrames = max(imageStruct.meta(3, :));
    imageFrame = cell(numFrames, 1);
    for i = 1:numFrames
        mask = imageStruct.meta(3, :) == i;
        imageFrame{i} = imageStruct.image(306:638, 334:717, mask); % WARNING: FXIME
    end
    
    % Infer from i-th frame to i+1-th frame
    quasiPointsHist = cell(numFrames, 1);
    quasiPointsHist{1} = points;
    flowPointsHist = cell(numFrames, 1);
    for i = 1:8
        
        fprintf('Calculating from frame %d to frame %d...\n', i, i+1);
        % optical flow
        flowPoints = runOpticalFlow(imageFrame{i}, imageFrame{i+1}, quasiPointsHist{i});
        % quasi-newton again
        quasiPoints = RefineLocation(imageFrame{i+1}, flowPoints);
        
        % record data
        flowPointsHist{i+1} = flowPoints;
        quasiPointsHist{i+1} = quasiPoints;
    end
    
    % visualzie
    if visualize
        figure
        hold on
        for i = 1:numFrames-1
            pts = quasiPointsHist{i};
            scatter3(pts(:, 1), pts(:, 2), pts(:, 3));
        end
    end
end
