function [ output_args ] = hd_stability()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
nSpikes=sum(data_video_spk(:,6));

% 6 degree bins
da=pi/30;
angBins=da/2:da:2*pi-da/2;
% Occupancy
histAng=hist(data_video_nospk(:,4),angBins);
%creating variable to check if all directional bins were sampled.
if length(find(histAng))==60
    bins_Sampled=1;
else
    bins_Sampled=0;
end
% Number of spikes per bin
spkPerAng=hist(data_video_spk(data_video_spk(:,6)==1,4),angBins);

% Tuning
hdTuning=(spkPerAng./histAng)*sampleRate;
% remove nan & inf
hdTuning(isnan(hdTuning) | isinf(hdTuning))=0;
data.(ratID{end-1}).(['S',strjoin(regexp(ratID{end},'\d*','Match'),'')]).hdTuning{i,event}=hdTuning;


end

