load('variables');

% [xHist, fvalHist, flagHist, debugHist] = GlobalImage(pt1(1:25, 1:30, :)); % run pt2 (one dot)
% [xHist, fvalHist, flagHist, debugHist] = GlobalImage(pt1); % run pt1
% CrunchyResult = GlobalImage(pt1);  % run pt3 (whole image)

pt3 = p(306:638, 334:717, :);
thred = quantile(pt3(:), 0.97)
pt3(pt3 > thred) = thred;
CrunchyResult = GlobalImage(pt3);
save('CrunchyResult.mat', 'CrunchyResult');
