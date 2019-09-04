function [latency, lat_coords]= time2Goal(x_smooth,y_smooth,rIN,rOut,mkplot,framerate,targetQuad)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%input:
%       -x_smooth: smooth x-coordinates (column vector)
%       -y_smooth: smooth y-coordinates (column vector)
%       -targetQuad: location of target of target quadrant in string form.
%           i.e. 'x,y' indicates target quadrant is within the positive
%           domains of x and y.
%       -frameRate: scalar in frames/s
%       -rOut: radius of pool (cm). 
%       -rIN: radius of inner cannulus (cm)
%       -mkplot: string argument. 'on' will plot coordinates outside of
%       annulus. 'off' does nothing. 
%output:
%       -time_quadT: Total 
%       -time_quadT: Total 
%       -time_quadT: Total 
%       -time_quadT: Total 
Lat_x=[];
Lat_y=[];
%Create function for pool diameter for plots
cosFxnOut = cos(0:2*pi/1000:2*pi)*rOut;
sinFxnOut = sin(0:2*pi/1000:2*pi)*rOut;

%create function for inner annulus
cosFxnIN = cos(0:2*pi/1000:2*pi)*rIN;
sinFxnIN = sin(0:2*pi/1000:2*pi)*rIN;

target= [cosFxnOut' sinFxnOut'; flip(cosFxnIN') flip(sinFxnIN')];

if strcmp(targetQuad,'x,y');
    target(target(:,1)<=0,:)=[]; target(target(:,2)<=0,:)=[];
elseif strcmp(targetQuad,'-x,-y');
    target(target(:,1)>=0,:)=[]; target(target(:,2)>=0,:)=[];
elseif strcmp(targetQuad,'-x,y');
    target(target(:,1)>=0,:)=[]; target(target(:,2)<=0,:)=[];
elseif strcmp(targetQuad,'x,-y');
    target(target(:,1)<=0,:)=[]; target(target(:,2)>=0,:)=[];
end
%create index of points inside inner annulus. 
[in,~]=inpolygon(x_smooth,y_smooth,target(:,1),target(:,2));

 i=1;
 while ismember(x_smooth(i,1),x_smooth(~in)) && ismember(y_smooth(i,1),y_smooth(~in))
     Lat_x=[Lat_x; x_smooth(i,1)];
     Lat_y=[Lat_y; y_smooth(i,1)];
     i=i+1;
 end
 lat_coords=[Lat_x Lat_y];
 latency=size(Lat_x,1)/framerate;
%make plots
if strcmp(mkplot,'on');
    hold on;
    title('Path to Goal')
    plot(cosFxnOut,sinFxnOut); hold on; %plot outside pool boundary
    plot(cosFxnIN,sinFxnIN);hold on; %plot inside cannulus 
    plot(x_smooth,y_smooth); hold on; %plot complete path
    plot(Lat_x,Lat_y,'*r'); %plot data points used in the lat calculation
    %plot quadrant boundaries
    n=4;
    tet=linspace(-pi,pi,n+1);
    xi=rOut*cos(tet)+0;
    yi=rOut*sin(tet)+0;
    for k=1:numel(xi)
        plot([0 xi(k)],[0 yi(k)]) %creates xy boundaries
        hold on
    end
end

