function [d,c]=DisplacementXY(m1,m2,s1,s2)
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

animation=1;

[RateMap1,~,~,~,~] = bindata(m1,30,s1,'no',76.5); RateMap1(isnan(RateMap1))=0; RateMap1(isinf(RateMap1))=0;
rALL=zeros(360,1);

if animation==1; figure(1); subplot(1,3,1);  pcolor(RateMap1); colormap jet; shading interp; end
for theta=0:359
    % ROTATE PATH & SPIKE XY
    XY=[cosd(theta) -sind(theta); sind(theta) cosd(theta)]*[m2(:,2),m2(:,3)]';
    spikes=[cosd(theta) -sind(theta); sind(theta) cosd(theta)]*[s2(:,2),s2(:,3)]';
    % CREATE RATEMAP
    try
    [RateMap,~,~,~,~] = bindata([[1:length(XY)]',XY'],30,[[1:size(spikes,2)]',spikes'],'no',76.5); RateMap(isnan(RateMap))=0;RateMap(isinf(RateMap))=0;
    catch
        test=1
    end
    % CORRELATE STANDARD AND ROTATED MAPS
    rALL(theta+1,1) = corr2(RateMap1,RateMap);
    if animation==1;subplot(1,3,2); pcolor(RateMap); colormap jet; shading interp;subplot(1,3,3); plot(rALL);pause(.25);end
end
% FIND MAX CORRELATION & DEGREE SHIFT
[c,d]=max(rALL);
end