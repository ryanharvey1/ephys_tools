function [gridout] = GridScore(SpatialAutoCorr,nBins)
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
%       - rGridScore
%       - maxSinuGrid
%       - meanSinuGrid
%       - SinuGrid
%       - rALL
%       - centerCUT
% 
% Created by Ben C Feb 2014
% Modified by Shawn W June 2014

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
rSinuScores = [];
GS1 = [];
SinuGrid = [];

for i = 1:length(list);
    centerCUT = SpatialAutoCorr;
    [rr cc] = meshgrid(1:(nBins*2)-1);
    circBINmid = sqrt((rr-nBins).^2+(cc-nBins).^2)<=((round(Minima/2)));
    circBINout = sqrt((rr-nBins).^2+(cc-nBins).^2)<=((round(Minima/2))+i);
    circBINout = ~circBINout;
    centerCUT(circBINout) = NaN;
    centerCUT(circBINmid) = NaN;

    % rotate by 0 degrees
    ai0 = imrotate(centerCUT,0,'nearest','crop');
    cCORR0 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR0 = reshape(ai0,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR0i = find(isnan(cCORR0)); % remove bins with zeros
    aiCORR0(cCORR0i,:) = [];
    cCORR0(cCORR0i,:) = [];
    aiCORR0i = find(isnan(aiCORR0));
    aiCORR0(aiCORR0i,:) = [];
    cCORR0(aiCORR0i,:) = [];
    r0 = corrcoef(cCORR0,aiCORR0); % obtain pearsons correlation between original and 6 deg rotated 'donut'
    
    % rotate by 6 degrees
    ai6 = imrotate(centerCUT,6,'nearest','crop');
    cCORR6 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR6 = reshape(ai6,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR6i = find(isnan(cCORR6)); % remove bins with zeros
    aiCORR6(cCORR6i,:) = [];
    cCORR6(cCORR6i,:) = [];
    aiCORR6i = find(isnan(aiCORR6));
    aiCORR6(aiCORR6i,:) = [];
    cCORR6(aiCORR6i,:) = [];
    r6 = corrcoef(cCORR6,aiCORR6); % obtain pearsons correlation between original and 6 deg rotated 'donut'
    
    % rotate by 12 degrees
    ai12 = imrotate(centerCUT,12,'nearest','crop');
    cCORR12 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR12 = reshape(ai12,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR12i = find(isnan(cCORR12)); % remove bins with zeros
    aiCORR12(cCORR12i,:) = [];
    cCORR12(cCORR12i,:) = [];
    aiCORR12i = find(isnan(aiCORR12));
    aiCORR12(aiCORR12i,:) = [];
    cCORR12(aiCORR12i,:) = [];
    r12 = corrcoef(cCORR12,aiCORR12); % obtain pearsons correlation between original and 12 deg rotated 'donut'
    
    % rotate by 18 degrees
    ai18 = imrotate(centerCUT,18,'nearest','crop');
    cCORR18 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR18 = reshape(ai18,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR18i = find(isnan(cCORR18)); % remove bins with zeros
    aiCORR18(cCORR18i,:) = [];
    cCORR18(cCORR18i,:) = [];
    aiCORR18i = find(isnan(aiCORR18));
    aiCORR18(aiCORR18i,:) = [];
    cCORR18(aiCORR18i,:) = [];
    r18 = corrcoef(cCORR18,aiCORR18); % obtain pearsons correlation between original and 18 deg rotated 'donut'

    % rotate by 24 degrees
    ai24 = imrotate(centerCUT,24,'nearest','crop');
    cCORR24 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR24 = reshape(ai24,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR24i = find(isnan(cCORR24)); % remove bins with zeros
    aiCORR24(cCORR24i,:) = [];
    cCORR24(cCORR24i,:) = [];
    aiCORR24i = find(isnan(aiCORR24));
    aiCORR24(aiCORR24i,:) = [];
    cCORR24(aiCORR24i,:) = [];
    r24 = corrcoef(cCORR24,aiCORR24); % obtain pearsons correlation between original and 24 deg rotated 'donut'

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

    % rotate by 36 degrees
    ai36 = imrotate(centerCUT,36,'nearest','crop');
    cCORR36 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR36 = reshape(ai36,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR36i = find(isnan(cCORR36)); % remove bins with zeros
    aiCORR36(cCORR36i,:) = [];
    cCORR36(cCORR36i,:) = [];
    aiCORR36i = find(isnan(aiCORR36));
    aiCORR36(aiCORR36i,:) = [];
    cCORR36(aiCORR36i,:) = [];
    r36 = corrcoef(cCORR36,aiCORR36); % obtain pearsons correlation between original and 36 deg rotated 'donut'

    % rotate by 42 degrees
    ai42 = imrotate(centerCUT,42,'nearest','crop');
    cCORR42 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR42 = reshape(ai42,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR42i = find(isnan(cCORR42)); % remove bins with zeros
    aiCORR42(cCORR42i,:) = [];
    cCORR42(cCORR42i,:) = [];
    aiCORR42i = find(isnan(aiCORR42));
    aiCORR42(aiCORR42i,:) = [];
    cCORR42(aiCORR42i,:) = [];
    r42 = corrcoef(cCORR42,aiCORR42); % obtain pearsons correlation between original and 42 deg rotated 'donut'
    
    % rotate by 48 degrees
    ai48 = imrotate(centerCUT,48,'nearest','crop');
    cCORR48 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR48 = reshape(ai48,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR48i = find(isnan(cCORR48)); % remove bins with zeros
    aiCORR48(cCORR48i,:) = [];
    cCORR48(cCORR48i,:) = [];
    aiCORR48i = find(isnan(aiCORR48));
    aiCORR48(aiCORR48i,:) = [];
    cCORR48(aiCORR48i,:) = [];
    r48 = corrcoef(cCORR48,aiCORR48); % obtain pearsons correlation between original and 48 deg rotated 'donut'

    % rotate by 54 degrees
    ai54 = imrotate(centerCUT,54,'nearest','crop');
    cCORR54 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR54 = reshape(ai54,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR54i = find(isnan(cCORR54)); % remove bins with zeros
    aiCORR54(cCORR54i,:) = [];
    cCORR54(cCORR54i,:) = [];
    aiCORR54i = find(isnan(aiCORR54));
    aiCORR54(aiCORR54i,:) = [];
    cCORR54(aiCORR54i,:) = [];
    r54 = corrcoef(cCORR54,aiCORR54); % obtain pearsons correlation between original and 54 deg rotated 'donut'
    
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

    % rotate by 66 degrees
    ai66 = imrotate(centerCUT,66,'nearest','crop');
    cCORR66 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR66 = reshape(ai66,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR66i = find(isnan(cCORR66)); % remove bins with zeros
    aiCORR66(cCORR66i,:) = [];
    cCORR66(cCORR66i,:) = [];
    aiCORR66i = find(isnan(aiCORR66));
    aiCORR66(aiCORR66i,:) = [];
    cCORR66(aiCORR66i,:) = [];
    r66 = corrcoef(cCORR66,aiCORR66); % obtain pearsons correlation between original and 66 deg rotated 'donut'
    
    % rotate by 72 degrees
    ai72 = imrotate(centerCUT,72,'nearest','crop');
    cCORR72 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR72 = reshape(ai72,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR72i = find(isnan(cCORR72)); % remove bins with zeros
    aiCORR72(cCORR72i,:) = [];
    cCORR72(cCORR72i,:) = [];
    aiCORR72i = find(isnan(aiCORR72));
    aiCORR72(aiCORR72i,:) = [];
    cCORR72(aiCORR72i,:) = [];
    r72 = corrcoef(cCORR72,aiCORR72); % obtain pearsons correlation between original and 72 deg rotated 'donut'
    
    % rotate by 78 degrees
    ai78 = imrotate(centerCUT,78,'nearest','crop');
    cCORR78 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR78 = reshape(ai78,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR78i = find(isnan(cCORR78)); % remove bins with zeros
    aiCORR78(cCORR78i,:) = [];
    cCORR78(cCORR78i,:) = [];
    aiCORR78i = find(isnan(aiCORR78));
    aiCORR78(aiCORR78i,:) = [];
    cCORR78(aiCORR78i,:) = [];
    r78 = corrcoef(cCORR78,aiCORR78); % obtain pearsons correlation between original and 78 deg rotated 'donut'
    
    % rotate by 84 degrees
    ai84 = imrotate(centerCUT,84,'nearest','crop');
    cCORR84 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR84 = reshape(ai84,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR84i = find(isnan(cCORR84)); % remove bins with zeros
    aiCORR84(cCORR84i,:) = [];
    cCORR84(cCORR84i,:) = [];
    aiCORR84i = find(isnan(aiCORR84));
    aiCORR84(aiCORR84i,:) = [];
    cCORR84(aiCORR84i,:) = [];
    r84 = corrcoef(cCORR84,aiCORR84); % obtain pearsons correlation between original and 84 deg rotated 'donut'
    
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

    % rotate by 96 degrees
    ai96 = imrotate(centerCUT,96,'nearest','crop');
    cCORR96 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR96 = reshape(ai96,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR96i = find(isnan(cCORR96)); % remove bins with zeros
    aiCORR96(cCORR96i,:) = [];
    cCORR96(cCORR96i,:) = [];
    aiCORR96i = find(isnan(aiCORR96));
    aiCORR96(aiCORR96i,:) = [];
    cCORR96(aiCORR96i,:) = [];
    r96 = corrcoef(cCORR96,aiCORR96); % obtain pearsons correlation between original and 96 deg rotated 'donut'
    
    % rotate by 102 degrees
    ai102 = imrotate(centerCUT,102,'nearest','crop');
    cCORR102 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR102 = reshape(ai102,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR102i = find(isnan(cCORR102)); % remove bins with zeros
    aiCORR102(cCORR102i,:) = [];
    cCORR102(cCORR102i,:) = [];
    aiCORR102i = find(isnan(aiCORR102));
    aiCORR102(aiCORR102i,:) = [];
    cCORR102(aiCORR102i,:) = [];
    r102 = corrcoef(cCORR102,aiCORR102); % obtain pearsons correlation between original and 102 deg rotated 'donut'
    
    % rotate by 108 degrees
    ai108 = imrotate(centerCUT,108,'nearest','crop');
    cCORR108 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR108 = reshape(ai108,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR108i = find(isnan(cCORR108)); % remove bins with zeros
    aiCORR108(cCORR108i,:) = [];
    cCORR108(cCORR108i,:) = [];
    aiCORR108i = find(isnan(aiCORR108));
    aiCORR108(aiCORR108i,:) = [];
    cCORR108(aiCORR108i,:) = [];
    r108 = corrcoef(cCORR108,aiCORR108); % obtain pearsons correlation between original and 108 deg rotated 'donut'
    
    % rotate by 114 degrees
    ai114 = imrotate(centerCUT,114,'nearest','crop');
    cCORR114 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR114 = reshape(ai114,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR114i = find(isnan(cCORR114)); % remove bins with zeros
    aiCORR114(cCORR114i,:) = [];
    cCORR114(cCORR114i,:) = [];
    aiCORR114i = find(isnan(aiCORR114));
    aiCORR114(aiCORR114i,:) = [];
    cCORR114(aiCORR114i,:) = [];
    r114 = corrcoef(cCORR114,aiCORR114); % obtain pearsons correlation between original and 114 deg rotated 'donut'
    
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
    
    % rotate by 126 degrees
    ai126 = imrotate(centerCUT,126,'nearest','crop');
    cCORR126 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR126 = reshape(ai126,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR126i = find(isnan(cCORR126)); % remove bins with zeros
    aiCORR126(cCORR126i,:) = [];
    cCORR126(cCORR126i,:) = [];
    aiCORR126i = find(isnan(aiCORR126));
    aiCORR126(aiCORR126i,:) = [];
    cCORR126(aiCORR126i,:) = [];
    r126 = corrcoef(cCORR126,aiCORR126); % obtain pearsons correlation between original and 126 deg rotated 'donut'
    
    % rotate by 132 degrees
    ai132 = imrotate(centerCUT,132,'nearest','crop');
    cCORR132 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR132 = reshape(ai132,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR132i = find(isnan(cCORR132)); % remove bins with zeros
    aiCORR132(cCORR132i,:) = [];
    cCORR132(cCORR132i,:) = [];
    aiCORR132i = find(isnan(aiCORR132));
    aiCORR132(aiCORR132i,:) = [];
    cCORR132(aiCORR132i,:) = [];
    r132 = corrcoef(cCORR132,aiCORR132); % obtain pearsons correlation between original and 132 deg rotated 'donut'
    
    % rotate by 138 degrees
    ai138 = imrotate(centerCUT,138,'nearest','crop');
    cCORR138 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR138 = reshape(ai138,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR138i = find(isnan(cCORR138)); % remove bins with zeros
    aiCORR138(cCORR138i,:) = [];
    cCORR138(cCORR138i,:) = [];
    aiCORR138i = find(isnan(aiCORR138));
    aiCORR138(aiCORR138i,:) = [];
    cCORR138(aiCORR138i,:) = [];
    r138 = corrcoef(cCORR138,aiCORR138); % obtain pearsons correlation between original and 138 deg rotated 'donut'
    
    % rotate by 144 degrees
    ai144 = imrotate(centerCUT,144,'nearest','crop');
    cCORR144 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR144 = reshape(ai144,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR144i = find(isnan(cCORR144)); % remove bins with zeros
    aiCORR144(cCORR144i,:) = [];
    cCORR144(cCORR144i,:) = [];
    aiCORR144i = find(isnan(aiCORR144));
    aiCORR144(aiCORR144i,:) = [];
    cCORR144(aiCORR144i,:) = [];
    r144 = corrcoef(cCORR144,aiCORR144); % obtain pearsons correlation between original and 144 deg rotated 'donut'

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
    
    % rotate by 156 degrees
    ai156 = imrotate(centerCUT,156,'nearest','crop');
    cCORR156 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR156 = reshape(ai156,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR156i = find(isnan(cCORR156)); % remove bins with zeros
    aiCORR156(cCORR156i,:) = [];
    cCORR156(cCORR156i,:) = [];
    aiCORR156i = find(isnan(aiCORR156));
    aiCORR156(aiCORR156i,:) = [];
    cCORR156(aiCORR156i,:) = [];
    r156 = corrcoef(cCORR156,aiCORR156); % obtain pearsons correlation between original and 156 deg rotated 'donut'
    
    % rotate by 162 degrees
    ai162 = imrotate(centerCUT,162,'nearest','crop');
    cCORR162 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR162 = reshape(ai162,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR162i = find(isnan(cCORR162)); % remove bins with zeros
    aiCORR162(cCORR162i,:) = [];
    cCORR162(cCORR162i,:) = [];
    aiCORR162i = find(isnan(aiCORR162));
    aiCORR162(aiCORR162i,:) = [];
    cCORR162(aiCORR162i,:) = [];
    r162 = corrcoef(cCORR162,aiCORR162); % obtain pearsons correlation between original and 162 deg rotated 'donut'
    
    % rotate by 168 degrees
    ai168 = imrotate(centerCUT,168,'nearest','crop');
    cCORR168 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR168 = reshape(ai168,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR168i = find(isnan(cCORR168)); % remove bins with zeros
    aiCORR168(cCORR168i,:) = [];
    cCORR168(cCORR168i,:) = [];
    aiCORR168i = find(isnan(aiCORR168));
    aiCORR168(aiCORR168i,:) = [];
    cCORR168(aiCORR168i,:) = [];
    r168 = corrcoef(cCORR168,aiCORR168); % obtain pearsons correlation between original and 168 deg rotated 'donut'
    
    % rotate by 174 degrees
    ai174 = imrotate(centerCUT,174,'nearest','crop');
    cCORR174 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR174 = reshape(ai174,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR174i = find(isnan(cCORR174)); % remove bins with zeros
    aiCORR174(cCORR174i,:) = [];
    cCORR174(cCORR174i,:) = [];
    aiCORR174i = find(isnan(aiCORR174));
    aiCORR174(aiCORR174i,:) = [];
    cCORR174(aiCORR174i,:) = [];
    r174 = corrcoef(cCORR174,aiCORR174); % obtain pearsons correlation between original and 174 deg rotated 'donut'
    
    % rotate by 180 degrees
    ai180 = imrotate(centerCUT,180,'nearest','crop');
    cCORR180 = reshape(centerCUT,((nBins*2)-1)*((nBins*2)-1),1); % reshape data so it is a column vector
    aiCORR180 = reshape(ai180,((nBins*2)-1)*((nBins*2)-1),1);
    cCORR180i = find(isnan(cCORR180)); % remove bins with zeros
    aiCORR180(cCORR180i,:) = [];
    cCORR180(cCORR180i,:) = [];
    aiCORR180i = find(isnan(aiCORR180));
    aiCORR180(aiCORR180i,:) = [];
    cCORR180(aiCORR180i,:) = [];
    r180 = corrcoef(cCORR180,aiCORR180); % obtain pearsons correlation between original and 180 deg rotated 'donut'
    
    % row of correlations between rotated and original 'donuts'
    r30s = [r30(1,2) r60(1,2) r90(1,2) r120(1,2) r150(1,2)];
    rGridScore = [rGridScore; r30s];
    rALL = [r0(1,2) r6(1,2) r12(1,2) r18(1,2) r24(1,2) r30(1,2) r36(1,2) r42(1,2) r48(1,2) r54(1,2) r60(1,2) r66(1,2) r72(1,2) r78(1,2) r84(1,2) r90(1,2) r96(1,2) r102(1,2) r108(1,2) r114(1,2) r120(1,2) r126(1,2) r132(1,2) r138(1,2) r144(1,2) r150(1,2) r156(1,2) r162(1,2) r168(1,2) r174(1,2) r180(1,2)];
    rSinuScores = [rSinuScores; rALL];
    
    % calculate min and max value for for group1 rotations (60, 120 degrees)
    Group1_GS = [r60(1,2) r120(1,2)];
    minGroup1_GS = min(Group1_GS);
    
    % calculate min and max value for group2 rotations (30, 90, 150 degrees)
    Group2_GS = [r30(1,2) r90(1,2) r150(1,2)];
    maxGroup2_GS = max(Group2_GS);
    
    % calculate grid score as difference between lowest maximum at 60 and
    % 120 and highest minimum at 30, 90, and 150 as in Sargolini et al 2006
    GS1_calc = minGroup1_GS - maxGroup2_GS;
    GS1 = [GS1; GS1_calc];
    
    % calculate grid score as fit of a sinusoid to the series of 6 degree
    % donut rotations
    time=[1:(length(rALL))]';
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Lower = [-Inf -Inf -Inf 0.5999];
    opts.StartPoint = [0 0 0 0.6];
    opts.Upper = [Inf Inf Inf 0.6];
    [FO, G]=fit(time,rALL','fourier1',opts);
    SinuGridi=G.rsquare;
    SinuGrid = [SinuGrid; SinuGridi];

    % calculate grid score based on buetfering et al 2014 and Allen et al 2014
    GS2 = ((rGridScore(:,2)+rGridScore(:,4))/2) - ((rGridScore(:,1)+rGridScore(:,3)+rGridScore(:,5))/3);
end

maxGS1 = max(GS1');
maxGS2 = max(GS2');
maxSinuGrid = max(SinuGrid');
meanSinuGrid = mean(SinuGrid');

% ADD TO STRUCT FOR OUTPUT
gridout.maxGS1=maxGS1;
gridout.maxGS2=maxGS2;
gridout.GS1=GS1;
gridout.GS2=GS2;
gridout.rGridScore=rGridScore;
gridout.maxSinuGrid=maxSinuGrid;
gridout.meanSinuGrid=meanSinuGrid;
gridout.SinuGrid=SinuGrid;
gridout.rALL=rALL;
gridout.centerCUT=centerCUT;