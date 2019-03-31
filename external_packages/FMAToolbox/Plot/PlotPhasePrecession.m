%PlotPhasePrecession - Plot phase precession plots and maps.
%
%  USAGE
%
%    PlotPhasePrecession(data,<options>,<options2>)
%
%    data           phase precession data obtained using <a href="matlab:help PhasePrecession">PhasePrecession</a>
%    <options>      optional list of property-value pairs (see table below)
%    <options2>     options for function <a href="matlab:help plot">plot</a>
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'      smoothing size in bins (0 = no smoothing, default = 4)
%     'nBins'       number of x and phase bins (default = [100 100])
%     'track'       'linear' or 'circular' (default = 'linear')
%    =========================================================================
%
%  SEE
%
%    See also PhasePrecession.

% Copyright (C) 2009-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function PlotPhasePrecession(data,varargin)

% Default values
track = 'linear';
nBins = [100 100];
smooth = 4;

% Check number of parameters
if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotPhasePrecession">PlotPhasePrecession</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PlotPhasePrecession">PlotPhasePrecession</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'smooth',
			smooth = varargin{i+1};
			if ~isdscalar(smooth,'>0'),
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help PlotPhasePrecession">PlotPhasePrecession</a>'' for details).');
			end
		case 'nbins',
			nBins = varargin{i+1};
			if ~isivector(nBins,'>0','#2'),
				error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help PlotPhasePrecession">PlotPhasePrecession</a>'' for details).');
			end
		case 'track',
			track = lower(varargin{i+1});
			if ~isstring(track,'linear','circular'),
				error('Incorrect value for property ''track'' (type ''help <a href="matlab:help PlotPhasePrecession">PlotPhasePrecession</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help PlotPhasePrecession">PlotPhasePrecession</a>'' for details).']);
	end
end

% Initialize a few variables
t = data.position.t;
x = data.position.x;
pp = data.position.phase;
r = data.rate.r;
pr = data.rate.phase;
nBinsX = nBins(1);
nBinsP = nBins(2);

% Setup axes
initialPosition = get(gca,'position');
delete(gca);
dy = initialPosition(4)/3;
position = initialPosition;
position(4) = initialPosition(4)/4;

% Phase vs position map
a = axes('position',position);
if strcmp(track,'linear'),
	plot([x;x],[pp;pp+2*pi]*180/pi,'k.','markersize',8);
else
	plot([x;x;x+1;x+1],[pp;pp+2*pi;pp;pp+2*pi]*180/pi,'k.','markersize',8);
end
set(a,'ytick',0:90:720,'ylim',[0 720],'tickdir','out');

% Phase vs position plot
position(2) = position(2)+dy;
a = axes('position',position);
x0 = Bin(x,[0 1],nBinsX);
p0 = Bin(pp,[0 2*pi],nBinsP);
map = Accumulate([p0 x0]);
if strcmp(track,'linear'),
	map = [map;map];
	map = Smooth(map,2);
	PlotColorMap(map,1,'x',linspace(0,1,nBinsX+1),'y',0:6:720);
else
	map = [map map;map map];
	map = Smooth(map,smooth);
	PlotColorMap(map,1,'x',linspace(0,2,2*nBinsX+1),'y',0:6:720);
end
set(a,'ytick',0:90:720,'ylim',[0 720],'tickdir','out');

% Phase vs rate plot
position(2) = position(2)+dy;
a = axes('position',position);
if isempty(r), return; end
rnd = randn(size(r))/5; % jitter points horizontally
plot(r+rnd,pr*180/pi,'k.','markersize',8);
hold on;
n = max(r);
for i = 1:n,
	if any(r==i),
		[m(i,1),c(i,:)] = CircularConfidenceIntervals(pr(r==i));
	else
		m(i,1) = NaN;
		c(i,1:2) = NaN;
	end
end
neg = m<0;
m(neg) = m(neg)+2*pi;
c(neg,:) = c(neg,:)+2*pi;
m = m*180/pi;
c = c*180/pi;
errorbar(1:n,m,m-c(:,2),c(:,1)-m,'r');
set(a,'ytick',0:90:360,'ylim',[0 360],'tickdir','out');
