%PlotSpikeWaveforms - Plot spike waveforms.
%
%  USAGE
%
%    p = PlotSpikeWaveforms(W,<options>,<options2>)
%
%    W              waveforms obtained using <a href="matlab:help GetSpikeWaveforms">GetSpikeWaveforms</a>
%    <options>      optional list of property-value pairs (see table below)
%    <options2>     options for function <a href="matlab:help plot">plot</a>
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'spacing'     vertical spacing between variables (default = 1 for point
%                   processes, 0 for continuous data)
%    =========================================================================
%
%  SEE
%
%    See also GetSpikeWaveforms, LoadSpikeWaveforms.

% Copyright (C) 2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function p = PlotSpikeWaveforms(W,varargin)

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help PlotSpikeWaveforms">PlotSpikeWaveforms</a>'' for details).');
end

% Default values
spacing = [];
v = {};

% Parse parameter list
n = 1;
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help PlotSpikeWaveforms">PlotSpikeWaveforms</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'spacing',
			spacing = varargin{i+1};
			if ~isdscalar(spacing),
				error('Incorrect value for property ''spacing'' (type ''help <a href="matlab:help PlotSpikeWaveforms">PlotSpikeWaveforms</a>'' for details).');
			end
		otherwise,
			v = varargin{i:end};
			if ~isa(v,'cell'), v = {v}; end
			break;
	end
end

if isempty(spacing),
	spacing = 50000;
end

% Plot
hold on;
n = size(W,2);
for i = 1:n,
	plot(spacing*(i-1)+permute(squeeze(W(:,i,:)),[2 1]),v{:});
end
