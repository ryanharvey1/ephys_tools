function [upsam, fig]=PerfectCircRateMap(ratemap,Plot)
% PerfectCircRateMap plot circular ratemap as a perfect circle
% 
%   Input:      
%           ratemap: n x p binned ratemap 
%           Plot: 1 or 0
% 
%   Notes:
%           This function requires the image processing toolbox for 'imresize'
% 
% Ryan E Harvey 2017

% DIM OF OUTPUT RATEMAP
x = -500:500;
y = -500:500;
[xx,yy] = meshgrid(x,y);
u = zeros(size(xx));
u((xx.^2+yy.^2)<500^2)=1;   % radius 500, center at the origin

% weight the points: point itself; average of nearest neighbors;
% averaged of diagonal neighbors.  These must add up to 1.
wp = .4;  wn = .4;  wd = .2;
ind = 2:length(x)-1;
u(ind,ind) = wp*u(ind,ind) ...
  + (wn/4)*(u(ind-1,ind  ) + u(ind+1,ind  ) + u(ind  ,ind-1) + u(ind  ,ind+1) ) ...
  + (wd/4)*(u(ind-1,ind-1) + u(ind-1,ind+1) + u(ind+1,ind-1) + u(ind+1,ind+1) );

% UPSAMPLE RATEMAP 
upsam=imresize(ratemap,[1001,1001]);

% CONVERT OUTSIDE CIRCLE BINS TO NAN
upsam(u<=0)=NaN;

% PLOT FIGURE
if Plot==1
    fig=figure; fig.Color=[1 1 1];
    h=pcolor(upsam);
    box off
    axis image
    axis off
    shading interp
end
end
