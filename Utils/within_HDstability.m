function [within_Coeff,within,normWithin] = within_HDstability(data_video_spk,sampleRate)
%Within_HDstability computes the 4-quarter stability score for head direction signals based off
%Boccara et al. 

%INPUT: 
%       -data_video_spk: timestamps interpolated with spike timestamps
%       -data_video_nospk: initial timestamps for occupancy 
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

time=[.25,.5,.75,1];

block=[data_video_spk(1,1),data_video_spk(end,1)*time];

for i=1:length(block)-1
    tempSpk=data_video_spk(data_video_spk(:,1)>block(i) & data_video_spk(:,1)<block(i+1),:);
    
    [~,~,~,~,~,hdTuning]=tuningcurve(tempSpk(tempSpk(:,6)==0,4),tempSpk(tempSpk(:,6)==1,4),sampleRate);

    within.hdTuning{i,1}=hdTuning;
        
%     figure;
%     plot(tempSpk(:,2),tempSpk(:,3),'.k');hold on
%     scatter(tempSpk(tempSpk(:,6)==1,2),tempSpk(tempSpk(:,6)==1,3),'Filled','r');
end

first=corr2(within.hdTuning{1,1},within.hdTuning{2,1});
second=corr2(within.hdTuning{1,1},within.hdTuning{3,1});
third=corr2(within.hdTuning{1,1},within.hdTuning{4,1});
fourth=corr2(within.hdTuning{2,1},within.hdTuning{3,1});
fifth=corr2(within.hdTuning{2,1},within.hdTuning{4,1});
sixth=corr2(within.hdTuning{3,1},within.hdTuning{4,1});

if sum(isnan([first,second,third,fourth,fifth,sixth]))>1
    test=1;
end

normTemp=[];
for i=1:4
    tempTune=within.hdTuning{i,:};
    normTemp=[normTemp; rescale(tempTune,0,1)];
end
normWithin=normTemp;
within_Coeff=nanmean([first,second,third,fourth,fifth,sixth]);
end

