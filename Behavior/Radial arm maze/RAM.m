% function RAM(X,Y)
%
% end
% RAM dims
pixelDist = 100/(range(data_video_smoothfilt2(:,2)));
ArmL=38/pixelDist;
CenterDi=24/pixelDist;

figure(1); plot(data_video_smoothfilt2(:,2),data_video_smoothfilt2(:,3))
hold on

xmax=max(data_video_smoothfilt2(:,2));
xmin=min(data_video_smoothfilt2(:,2));
ymax=max(data_video_smoothfilt2(:,3));
ymin=min(data_video_smoothfilt2(:,3));

center=[xmax-xmin,ymax-ymin];

 x=median([xmax,xmin]); y=median([ymax,ymin]); rad=CenterDi/2;
th = 0:pi/179.5:2*pi; % 0 to 2*pi(6.28318530717959) at 0.0175 increments to equal 360 points
xunit = rad * cos(th) + x; yunit = rad * sin(th) + y;


figure(1); plot(xunit,yunit)