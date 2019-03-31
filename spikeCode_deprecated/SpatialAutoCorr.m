function [SpatialAutoCorr] = SpatialAutoCorr(v,nBins)
% SpatialAutoCorr generates a spatial autocorrelogram for a rate map
%
% Input: 
%       - v = array of firing rates
%
% Output:
%       - SpatialAutoCorr = array of correlation coefficients
% 
% Ben C 2014

% spatial autocorr computed using function 
% from http://www.mathworks.com/matlabcentral/fileexchange/29005-generalized-normalized-cross-correlation 
v(isnan(v))=0; v(isinf(v))=0;
[c,numberOfOverlapPixels] = normxcorr2_general(v,v);

% reshape autocorr and #ofpixels arrays so they are a column vector
ci = reshape(c,((nBins*2)-1)*((nBins*2)-1),1);
pixelsi = reshape(numberOfOverlapPixels,((nBins*2)-1)*((nBins*2)-1),1);

% keep only bins with >20 samples
Twentyci = find(pixelsi < 20);
ci(Twentyci,:) = NaN;
SpatialAutoCorr = reshape(ci,(nBins*2)-1,(nBins*2)-1);