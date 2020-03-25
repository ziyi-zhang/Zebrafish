function [sigmaRes] = SubtractionSigma(xyArray, method)
% [sigmaRes] = SubtractionSigma(xyArray[, method])
% Returns the subraction function sigma based on input sample points
% 'xyArray' 2*N array of sample location

    if nargin<2
        method = 'hard';
    end
    
    if strcmp(method, 'hard')
        filter = (xyArray(1, :).^2 + xyArray(2, :).^2) > 1;
        sigmaRes = ones(1, size(xyArray, 2));
        sigmaRes(filter) = -1;
    end
    
    if strcmp(method, 'sigmoid')
        dist = sqrt(xyArray(1, :).^2 + xyArray(2, :).^2);
        sigmaRes = 2 ./ (1 + exp(6.*dist-6)) - 1;
    end
end
