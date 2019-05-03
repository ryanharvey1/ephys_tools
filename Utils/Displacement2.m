function [d,c]=Displacement2(RateMap1,RateMap2)
% DisplacementXY
%   Input:
%           m1: xy from first path [t,x,y]
%           m2: xy from second [t,x,y]
%           s1: xy from first spikes [t,x,y]
%           s2: xy from second spikes [t,x,y]
%   Output:
%           d: Degree of rotation in field
%           c: Correlation to sample ratemap given rotation
%
%   other:
%           setting animation to 1 will plot the rotation and correlations
% Ryan E. Harvey 2017

animation=0;

% [RateMap1,~,~,~,~] = bindata(m1,30,s1,'no',76.5);
RateMap1(isnan(RateMap1))=0; RateMap1(isinf(RateMap1))=0;
RateMap1 = padarray(RateMap1,[3 3],0,'both');

RateMap2(isnan(RateMap2))=0; RateMap2(isinf(RateMap2))=0;
RateMap2 = padarray(RateMap2,[3 3],0,'both');

rALL=zeros(360,1);

if animation==1; figure(1); subplot(1,3,1);  imagesc(RateMap1); colormap jet; end
for theta=0:359
    % ROTATE PATH & SPIKE XY
    RateMap=imrotate(RateMap2,theta,'nearest','crop');

    % CORRELATE STANDARD AND ROTATED MAPS
    rALL(theta+1,1) = corr2(RateMap1,RateMap);
    if animation==1;subplot(1,3,2); imagesc(RateMap); colormap jet; subplot(1,3,3); plot(rALL);pause(.25);end
end
% FIND MAX CORRELATION & DEGREE SHIFT
[c,d]=max(rALL);
end