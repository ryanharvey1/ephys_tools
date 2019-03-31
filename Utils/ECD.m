function [pval,r]=ECD(tsxy,headang,spkts,cuexy)
% ECD ego-centric cue direction (Wilber et al. 2014)
%
% WORKING
addpath /Users/ryanharvey/GoogleDrive/MatlabDir/CircStat2012a

if isempty(headang)
    headang=XYangle(tsxy(:,2),tsxy(:,3));
    headang=[headang;headang(end)];
%     ang=wrapTo180(ang);
end

ecdangle=atan2d(cuexy(2)-tsxy(:,3),cuexy(1)-tsxy(:,2));

angle=wrapTo180(headang-ecdangle);


% CREATE TUNING CURVE
% 6 degree bins
angBins=linspace(-pi,pi,60);
% Occupancy
histAng=hist(deg2rad(angle),angBins);
% Number of spikes per bin
spkPerAng=hist(deg2rad(angle(logical(tsxy(:,end)))),angBins);
% Tuning
hdTuning=(spkPerAng./histAng)*30;
% remove nan & inf
hdTuning(isnan(hdTuning) | isinf(hdTuning))=0;

[pval,~] = circ_rtest(angBins',hdTuning',deg2rad(6));


% STABILITY
half=round(length(tsxy)/2);

histAng=hist(deg2rad(angle(1:half)),angBins);
spkPerAng=hist(deg2rad(angle(logical(tsxy(1:half,end)))),angBins);
hdTuning=(spkPerAng./histAng)*30;
hdTuning(isnan(hdTuning) | isinf(hdTuning))=0;
hdTuning1=hdTuning;

histAng=hist(deg2rad(angle(half+1:end)),angBins);
spkPerAng=hist(deg2rad(angle(logical(tsxy(half+1:end,end)))),angBins);
hdTuning=(spkPerAng./histAng)*30;
hdTuning(isnan(hdTuning) | isinf(hdTuning))=0;
hdTuning2=hdTuning;

r=corr2(hdTuning1,hdTuning2);
end

