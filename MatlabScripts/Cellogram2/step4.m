% output
% [xyzreHist_processed] = step4(xyzreHist, cropData);

function [xyzreHist_processed] = step4(xyzreHist, cropData)

     frames = length(xyzreHist);
     
     for f = 1:frames
        xyzreHist{f}(:, 1) = xyzreHist{f}(:, 1) + cropData.c0 - 1;  % x
        xyzreHist{f}(:, 2) = xyzreHist{f}(:, 2) + cropData.r0 - 1;  % y
        xyzreHist{f}(:, 3) = xyzreHist{f}(:, 3) + (cropData.sliceBegin - 1) + 1;  % bottom center with height 2
     end
     xyzreHist_processed = xyzreHist;
end
