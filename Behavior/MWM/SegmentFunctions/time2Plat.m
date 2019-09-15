function [latency, lat_coords,dwell_Time]= time2Plat(x_smooth,y_smooth,x_plat,y_plat,rPlat,rPool,mkplot,framerate)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%input:
%       -x_smooth: smooth x-coordinates (column vector)
%       -y_smooth: smooth y-coordinates (column vector)
%       -x_plat: x-coordinate of platform origin
%       -y_plat: y-coordinate of platform origin
%       -frameRate: scalar in frames/s
%       -rOut: radius of pool (cm). 
%       -rPlat: radius of platform (cm)
%       -mkplot: string argument. 'on' will plot coordinates outside of
%       annulus. 'off' does nothing. 
%output:
%       -latency: time to platform area
%       -lat_coords: latency coordinates
%       -dewll_Time: amount of time in seconds the rat was in the platform
%       area. 

Lat_x=[];
Lat_y=[];
%Create function for pool diameter for plots
cosFxnOut = cos(0:2*pi/1000:2*pi)*rPool;
sinFxnOut = sin(0:2*pi/1000:2*pi)*rPool;

%create function for platform boundaries
cosFxnIN = x_plat+cos(0:2*pi/1000:2*pi)*rPlat;
sinFxnIN = y_plat+sin(0:2*pi/1000:2*pi)*rPlat;

%create index of points inside inner annulus. 
[in,~]=inpolygon(x_smooth,y_smooth,cosFxnIN,sinFxnIN);

 i=1;
 while ismember(x_smooth(i,1),x_smooth(~in)) && ismember(y_smooth(i,1),y_smooth(~in))
     Lat_x=[Lat_x; x_smooth(i,1)];
     Lat_y=[Lat_y; y_smooth(i,1)];
     i=i+1;
     if i>size(x_smooth,1)
         break
     end
 end
 lat_coords=[Lat_x Lat_y];
 latency=size(Lat_x,1)/framerate;
 dwell_Time=length(x_smooth(in) & y_smooth(in))/framerate;
%make plots
figure;
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
    xi=rPool*cos(tet)+0;
    yi=rPool*sin(tet)+0;
    for k=1:numel(xi)
        plot([0 xi(k)],[0 yi(k)]) %creates xy boundaries
        hold on
    end
end
