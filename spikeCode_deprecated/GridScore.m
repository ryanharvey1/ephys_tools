function gridout = GridScore(SpatialAutoCorr,nBins)
% GridScore creates autocorr map and calcs grid score based on Stensola et al 2012
%
% Input: 
%       - SpatialAutoCorr = array of correlation coefficients
%
% Output:
%       - maxGS1
%       - maxGS2
%       - GS1
%       - GS2
%       - centerCUT
%       - rGridScore
% 
% Created by Ben C Feb 2014

% take mean of correlations surrounding center of autocorr
cCENTER = SpatialAutoCorr(nBins,:)';
cCENTER(:,2) = SpatialAutoCorr(:,nBins);
MeanCenter = mean(cCENTER');

% invert mean around center and find first peak (local minima around center of SpatialAutoCorr)
mCenterInv = 1.01 * max(MeanCenter') - MeanCenter';
[MinIdx] = findpeaks(mCenterInv);
Minima = min(abs((MinIdx(1).loc)-nBins));

% define ring around center circle and increase the size of the circle for
% each loop
StartSuccSamples = ((round(Minima/2)))+2;
maxSuccSamples = size(SpatialAutoCorr);
maxSuccSamples = maxSuccSamples(1,1) - 2;
list = StartSuccSamples:1:maxSuccSamples;
rGridScore = [];
GS1 = [];
GS2 = [];

for i = 1:length(list);
    centerCUT = SpatialAutoCorr;
    [rr cc] = meshgrid(1:(nBins*2)-1);
    circBINmid = sqrt((rr-nBins).^2+(cc-nBins).^2)<=((round(Minima/2)));
    circBINout = sqrt((rr-nBins).^2+(cc-nBins).^2)<=((round(Minima/2))+i);
    circBINout = ~circBINout;
    centerCUT(circBINout) = NaN;
    centerCUT(circBINmid) = NaN;
    
    % rotate by 30 degrees
    ai30 = imrotate(centerCUT,30,'nearest','crop');
    cCORR30 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR30 = reshape(ai30,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR30i = find(isnan(cCORR30)); % remove bins with zeros
    aiCORR30(cCORR30i,:) = [];
    cCORR30(cCORR30i,:) = [];
    aiCORR30i = find(isnan(aiCORR30));
    aiCORR30(aiCORR30i,:) = [];
    cCORR30(aiCORR30i,:) = [];
    r30 = corrcoef(cCORR30,aiCORR30); % obtain pearsons correlation between original and 30 deg rotated 'donut'

    % rotate by 60 degrees
    ai60 = imrotate(centerCUT,60,'nearest','crop');
    cCORR60 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR60 = reshape(ai60,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR60i = find(isnan(cCORR60)); % remove bins with zeros
    aiCORR60(cCORR60i,:) = [];
    cCORR60(cCORR60i,:) = [];
    aiCORR60i = find(isnan(aiCORR60));
    aiCORR60(aiCORR60i,:) = [];
    cCORR60(aiCORR60i,:) = [];
    r60 = corrcoef(cCORR60,aiCORR60); % obtain pearsons correlation between original and 60 deg rotated 'donut'

    % rotate by 90 degrees
    ai90 = imrotate(centerCUT,90,'nearest','crop');
    cCORR90 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR90 = reshape(ai90,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR90i = find(isnan(cCORR90)); % remove bins with zeros
    aiCORR90i = find(isnan(aiCORR90));
    aiCORR90(cCORR90i,:) = [];
    cCORR90(cCORR90i,:) = [];
    aiCORR90i = find(isnan(aiCORR90));
    aiCORR90(aiCORR90i,:) = [];
    cCORR90(aiCORR90i,:) = [];
    r90 = corrcoef(cCORR90,aiCORR90); % obtain pearsons correlation between original and 90 deg rotated 'donut'

    % rotate by 120 degrees
    ai120 = imrotate(centerCUT,120,'nearest','crop');
    cCORR120 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR120 = reshape(ai120,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR120i = find(isnan(cCORR120)); % remove bins with zeros
    aiCORR120i = find(isnan(aiCORR120));
    aiCORR120(cCORR120i,:) = [];
    cCORR120(cCORR120i,:) = [];
    aiCORR120i = find(isnan(aiCORR120));
    aiCORR120(aiCORR120i,:) = [];
    cCORR120(aiCORR120i,:) = [];
    r120 = corrcoef(cCORR120,aiCORR120); % obtain pearsons correlation between original and 120 deg rotated 'donut'

    % rotate by 150 degrees
    ai150 = imrotate(centerCUT,150,'nearest','crop');
    cCORR150 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR150 = reshape(ai150,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR150i = find(isnan(cCORR150)); % remove bins with zeros
    aiCORR150i = find(isnan(aiCORR150));
    aiCORR150(cCORR150i,:) = [];
    cCORR150(cCORR150i,:) = [];
    aiCORR150i = find(isnan(aiCORR150));
    aiCORR150(aiCORR150i,:) = [];
    cCORR150(aiCORR150i,:) = [];
    r150 = corrcoef(cCORR150,aiCORR150); % obtain pearsons correlation between original and 150 deg rotated 'donut'
    
    % row of correlations between rotated and original 'donuts'
    rALL = [r30(1,2) r60(1,2) r90(1,2) r120(1,2) r150(1,2)];
    rGridScore = [rGridScore; rALL];
    
    % calculate min and max value for for group1 rotations (60, 120 degrees)
    Group1_GS = [r60(1,2) r120(1,2)];
    maxGroup1_GS = max(Group1_GS);
    minGroup1_GS = min(Group1_GS);
    
    % calculate min and max value for group2 rotations (30, 90, 150 degrees)
    Group2_GS = [r30(1,2) r90(1,2) r150(1,2)];
    maxGroup2_GS = max(Group2_GS);
    minGroup2_GS = min(Group2_GS);
    
    % calculate grid score as difference between lowest maximum at 60 and
    % 120 and highest minimum at 30, 90, and 150 as in Sargolini et al 2006
    GS1_calc = minGroup1_GS - maxGroup2_GS;
    GS1 = [GS1; GS1_calc];
    
    % calculate grid score as difference between highest maximum at 60 and 120 and lowest
    % minimum at 30, 90, and 150 
    GS2_calc = maxGroup1_GS - minGroup2_GS;
    GS2 = [GS2; GS2_calc];
end

maxGS1 = max(GS1');
maxGS2 = max(GS2');

gridout.maxGS1=maxGS1;
gridout.maxGS2=maxGS2;
gridout.GS1=GS1; 
gridout.GS2=GS2;
gridout.rGridScore=rGridScore;
gridout.centerCUT=centerCUT;
