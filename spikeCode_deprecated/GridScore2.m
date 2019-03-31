function [maxSinuGrid] = GridScore2(SpatialAutoCorr,nBins)
% GridScore creates autocorr map and calcs grid score based on Stensola et al 2012
%
% Input: 
%       - SpatialAutoCorr = array of correlation coefficients
%
% Output:
%       - maxGS
%       - Group1_GS
%       - Group2_GS
%       - centerCUT
%       - rGridScore
% 
% Created by Ben C Feb 2014

% set options for sinusoidal fit
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Lower = [-Inf -Inf -Inf 0.5999];
opts.StartPoint = [0 0 0 0.6];
opts.Upper = [Inf Inf Inf 0.6];

% take mean of correlations surrounding center of autocorr
cCENTER = SpatialAutoCorr(nBins,:)';
cCENTER(:,2) = SpatialAutoCorr(:,nBins);
MeanCenter = mean(cCENTER');

% invert mean around center and find first peak (local minima around center of SpatialAutoCorr)
mCenterInv = 1.01 * max(MeanCenter') - MeanCenter';
[MinIdx] = findpeaks(mCenterInv);
Minima = min(abs((MinIdx(1))-nBins));

% define ring around center circle and increase the size of the circle for
% each loop
StartSuccSamples = ((round(Minima/2)))+2;
maxSuccSamples = size(SpatialAutoCorr);
maxSuccSamples = maxSuccSamples(1,1) - 2;
list = StartSuccSamples:1:maxSuccSamples;
dirBins = 0:6:180;
rALL = [];
rGridScore = [];

for i = 1:length(list)
    centerCUT = SpatialAutoCorr;
    [rr cc] = meshgrid(1:(nBins*2)-1);
    circBINmid = sqrt((rr-nBins).^2+(cc-nBins).^2)<=((round(Minima/2)));
    circBINout = sqrt((rr-nBins).^2+(cc-nBins).^2)<=((round(Minima/2))+i);
    circBINout = ~circBINout;
    centerCUT(circBINout) = NaN;
    centerCUT(circBINmid) = NaN;
    for j = 1:length(dirBins)
        ai = imrotate(centerCUT,j,'nearest','crop');
        cCORR = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
        aiCORR = reshape(ai,((nBins*2)-1)*((nBins*2)-1),1);
        cCORRi = find(isnan(cCORR)); % remove bins with NaNs
        aiCORR(cCORRi,:) = [];
        cCORR(cCORRi,:) = [];
        aiCORRi = find(isnan(aiCORR)); % remove bins with NaNs
        aiCORR(aiCORRi,:) = [];
        cCORR(aiCORRi,:) = [];
        r = corrcoef(cCORR,aiCORR); % obtain pearsons correlation between original and rotated 'donut'
        rALL(j)=r(1,2);
%         rALL = [rALL; r];
    end
    
    % calculate grid score as fit of a sinusoid to the series of 6 degree
    % donut rotations
    [~, G]=fit([1:(length(rALL))]',rALL','fourier1',opts);
    SinuGrid(i,1)=G.rsquare;
    
end
maxSinuGrid = max(SinuGrid);

%     rALL = rALL';
%     rGridScore = [rGridScore; rALL]; 

