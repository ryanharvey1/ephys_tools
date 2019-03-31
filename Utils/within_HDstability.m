function [within_Coeff,within,normWithin] = within_HDstability(data_video_spk,data_video_nospk,sampleRate,Angle,Spike)
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
idx=1;
sessionLength=length(data_video_nospk(:,1));

for i=1:sessionLength*.25:sessionLength
    tempNoSpk=data_video_nospk(round(i):round(sessionLength*time(idx)),:);
    tempSpk=data_video_spk(round(i):round(sessionLength*time(idx)),:);
    % 6 degree bins
    da=pi/30;
    angBins=da/2:da:2*pi-da/2;
    % Occupancy
    histAng=hist(tempNoSpk(:,Angle),angBins);
    
    % Number of spikes per bin
    spkPerAng=hist(tempSpk(tempSpk(:,Spike)==1,Angle),angBins);
    
    % Tuning
    hdTuning=(spkPerAng./histAng)*sampleRate;
    
    % remove nan & inf
    hdTuning(isnan(hdTuning) | isinf(hdTuning))=0;
    
    within.hdTuning{idx,1}=hdTuning;
    
    idx=idx+1;
    
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

