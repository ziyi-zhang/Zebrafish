
prevPoints = [score.meanX, score.meanY, score.meanZ];
pt3F2 = imageC1F2.image(306:638, 334:717, :);
resPoints = runOpticalFlow(pt3F1, pt3F2, prevPoints);
ImshowTiff(pt3F2(:, :, 15)+pt3F2(:, :, 18)+pt3F2(:, :, 19));
hold on
viscircles([prevPoints(:, 1), prevPoints(:, 2)], score.meanR, 'LineWidth', 0.5, 'Color', [0 0.4470 0.7410]);
viscircles([resPoints(:, 1), resPoints(:, 2)], score.meanR, 'LineWidth', 0.5, 'Color', [0.8500 0.3250 0.0980]);
%scatter(prevPoints(:, 1), prevPoints(:, 2))
%scatter(resPoints(:, 1), resPoints(:, 2))
