function [out] = SpatialAutoCorr(v,nBins)

com=which('SpatialAutoCorr');
com=strsplit(com,filesep);
basedir=[com{1},filesep,'Users',filesep,com{3},filesep,'GoogleDrive',filesep,'MatlabDir'];
addpath([basedir,filesep,filesep,'BClarkToolbox',filesep,'Analyses',filesep,'spikeCode'],...
    [basedir,filesep,'BClarkToolbox',filesep, 'Analysis']);

v(isnan(v))=0; v(isinf(v))=0;
[c,numberOfOverlapPixels] = normxcorr2_general(v,v);

% reshape autocorr and #ofpixels arrays so they are a column vector
ci = reshape(c,((nBins*2)-1)*((nBins*2)-1),1);
pixelsi = reshape(numberOfOverlapPixels,((nBins*2)-1)*((nBins*2)-1),1);

% keep only bins with >20 samples
Twentyci = find(pixelsi < 20);
ci(Twentyci,:) = NaN;
out = reshape(ci,(nBins*2)-1,(nBins*2)-1);
end