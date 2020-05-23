function [] = LoGParameter(image, hsizeArr, sigmaArr, zoomIn)
% [] = LoGParameter(image[, hsizeArr][, sigmaArr][, zoomIn])
% This function is used to visualize the effect of LoG (Laplacian of
% Gaussian) filter size and sigma.
% 'image' is a 2D matrix of a grayscale image
% 'hsizeArr' is the candidate filter size/length array
% 'sigmaArr' is the candidate Gaussian sigma array
% 'zoomIn' is the ratio of zoom in and appies on all subplots
% Example Use:
% LoGParameter(image, [], [], 2.5)
% Ziyi. Feb, 2020.

   if nargin<4, zoomIn = 1;end
   if nargin<3 || isempty(sigmaArr), sigmaArr = 1:4;end
   if nargin<2 || isempty(hsizeArr), hsizeArr = 4:2:16;end

   % zoom in calculation
   xzoomIn(1) = ceil(size(image, 1) / 2 * (1 - 1.0/zoomIn));
   xzoomIn(2) = size(image, 1) - xzoomIn(1);
   yzoomIn(1) = ceil(size(image, 2) / 2 * (1 - 1.0/zoomIn));
   yzoomIn(2) = size(image, 2) - yzoomIn(1);
   if xzoomIn(1)<1, xzoomIn(1)=1;end
   if yzoomIn(1)<1, yzoomIn(1)=1;end
   % plot
   figure
   ha = tight_subplot(length(sigmaArr), length(hsizeArr));
   count = 0;
   for sigma = sigmaArr
      for hsize = hsizeArr 
       
           count = count + 1;
           axes(ha(count)); %#ok
           h = fspecial('log', hsize, sigma);
           image_ = imfilter(image, h);

           image_ = image_(xzoomIn(1):xzoomIn(2), yzoomIn(1):yzoomIn(2)); % zoom in
           imshow(image_, [min(image_, [], 'all'), max(image_, [], 'all')]);
           colormap bone
           titleStr = "l:" + hsize + " s:" + sigma;
           title(titleStr);
       end
   end
end
