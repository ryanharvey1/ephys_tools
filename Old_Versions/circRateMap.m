function [ circBinMap ] = circRateMap(FilledRateMatrix)
%Circ Rate map converts the bins outside of your circle to NaNs for cell
%recordings in a cylinder. By L.Berk March 2017
%   Input:
%       -FilledRateMatrix - consists of unsmoothed ratematrix
%   Output: 
%       -rate matrix in the form of a circle 
% 
    [r,c]=size(FilledRateMatrix);
    %create circular mask
    [W,H] = meshgrid(1:c,1:r);
    mask = logical(((W-c/2).^2 + (H-r/2).^2) < (max([r,c])/2)^2);
    %convert outside of circle to NaNs
    circBinMap=FilledRateMatrix;
    circBinMap(~logical(reshape(mask,[],length(FilledRateMatrix))))=NaN;
end

