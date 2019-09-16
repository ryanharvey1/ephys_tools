function [r,I,Ispk,peakrate,prefdirec,hdTuning]=tuningcurve(a,spk_a,Fs)
% tuningcurve: creates a directional tuning curve and calculates several
% circular distribution measures 
%
% Input:
%       a: a vector of head-direction (degrees) from video-tracking 
%       spk_a: head-direction (degrees) at times of spikes
%       Fs: video sampling frequency 
% Output:
%       r: mean resultant vector length
%       I: Information rate (bit/sec)
%       Ispk: Information per spike (bit/spk)
%       peakrate: peak firing rate
%       prefdirec: preferred direction 
%       hdTuning: Tuning curve 
% Adapted from Adrien Peyrache github.com/PeyracheLab
% by Ryan Harvey (2018)

%angular bins
% da=pi/30;angBins=da/2:da:2*pi-da/2;
angBins=0:6:360;
bin_centers=movmedian(angBins,2);
bin_centers(1)=[];
% Occupancy
histAng=histcounts(a,angBins);
% Number of spikes per bin
spkPerAng=histcounts(spk_a,angBins);
% Tuning
hdTuning=(spkPerAng./histAng)*Fs;
% remove nan & inf
hdTuning(isnan(hdTuning) | isinf(hdTuning))=0;
% calculate mean resultant vector length
r=circ_r(deg2rad(bin_centers)',hdTuning',deg2rad(6));
% peak rate
[peakrate,I]=max(hdTuning);
% preferred direction
prefdirec=angBins(I);
  

%% Calculate information per second and information per spike
% Following Skaggs et al., 1993
% by Adrien Peyrache

% length of recording
T = size(a,1)/Fs;
% Number of spikes
N = size(spk_a,1);
% Average firing rate
fr = N/T;
% probability of occupancy:
Px = histAng./sum(histAng);
logTerm = log2(hdTuning/fr);
% Correct for undefined values
logTerm(hdTuning==0) = 0;
% Little trick to express a sum as a dot product
I = hdTuning * (logTerm.*Px)' ;
% Divide by firing rate to obtain information per spike 
Ispk = I/fr;

end