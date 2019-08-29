function [within_Coeff,within,normWithin] = within_HDstability(data_video_spk,sampleRate)
%Within_HDstability computes the 4-quarter stability score for head direction signals based off
%Boccara et al. 

%INPUT: 
%       -data_video_spk: timestamps interpolated with spike timestamps
%       -sampleRate: video sample rate (in hz)
%       -Angle: 
%       -Spike: 

%OUTPUT: 
%       -within_Coeff: 4-qrt stability score
%       -structure of raw tuning curves for each qtr in ascending order
%       (e.g. row 1= qtr 1). 
%       -normWithin: matrix of normalize tuning curves (row 1= qtr 1, row
%       2=qtr 2 etc.)

% Created by LBerkowitz March 2018, updated by LB July 2018 

time = [.25,.5,.75,1];

start = min(data_video_spk(:,1));
stop = max(data_video_spk(:,1));
blocking_factor = stop-start; 

block = blocking_factor*time;

clear blocking_factor start stop time

block = [data_video_spk(1,1) data_video_spk(1,1)+block];

for i = 1:length(block)-1
    tempSpk = data_video_spk(data_video_spk(:,1) > block(i) & data_video_spk(:,1) < block(i+1),:);
    
    [~,~,~,~,~,hdTuning] = tuningcurve(tempSpk(tempSpk(:,6)==0,4),tempSpk(tempSpk(:,6)==1,4),sampleRate);

    within.hdTuning{i,1} = hdTuning;
 
end
 clear block
 
 %calculate correlation for all pairs
first=corr2(within.hdTuning{1,1},within.hdTuning{2,1});
second=corr2(within.hdTuning{1,1},within.hdTuning{3,1});
third=corr2(within.hdTuning{1,1},within.hdTuning{4,1});
fourth=corr2(within.hdTuning{2,1},within.hdTuning{3,1});
fifth=corr2(within.hdTuning{2,1},within.hdTuning{4,1});
sixth=corr2(within.hdTuning{3,1},within.hdTuning{4,1});

%warning if data is nan
if sum(isnan([first,second,third,fourth,fifth,sixth]))>1
    warning('NaNs within the data are preventing an accurate calculation of stability. Check your data')
end

normTemp=[];
for i=1:4
    tempTune=within.hdTuning{i,:};
    normTemp=[normTemp; rescale(tempTune,0,1)];
end
normWithin=normTemp;
within_Coeff=nanmean([first,second,third,fourth,fifth,sixth]);
end

