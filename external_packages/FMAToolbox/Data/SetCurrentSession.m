%SetCurrentSession - Load all data for a given recording session.
%
% Set current session files and read data from disk. Calling SetCurrentSession
% without parameters will display a file selection dialog.
%
%  USAGE
%
%    SetCurrentSession(filename)
%
%    filename            optional parameter file name; use 'same' to force reload
%
%  NOTE
%
%    If no parameter file name is specified, an interactive file selection
%    dialog is displayed.

% Copyright (C) 2004-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function SetCurrentSession(filename)

global DATA;

% Initialization
if isempty(DATA),
	format long g;
	DATA.session.basename = '';
	DATA.session.path = '';
	DATA.positions = [];
	DATA.spikes = [];
	DATA.events = [];
	% Default settings
	Settings;
end

if nargin == 0 || (strcmp(filename,'same') && isempty(DATA.session.basename)),
	% Interactive mode
	[filename,path] = uigetfile('*.xml','Please select a parameter file for; this session');
	if filename == 0,return; end
	filename = [path filename];
end

if strcmp(filename,'same'),
	% Force reload
	path = DATA.session.path;
	basename = DATA.session.basename;
else
%  	if isempty(findstr(filename,'/')),
%  		filename = [pwd '/' filename];
%  	end
	% Parse file name
	while filename(end) == '/', filename(end) = []; end
	[path,basename] = fileparts(filename);
	if isempty(path)|path(1)~='/', path = [pwd '/' path]; end
	path = strrep(path,'/./','/');
end

disp(['Loading session files for ' basename]);

% File already loaded?
if strcmp(basename,DATA.session.basename) & strcmp(path,DATA.session.path) & ~strcmp(filename,'same'),
	disp(['... session files already loaded, skipping - type SetCurrentSession(''same'') to force reload']);
	disp('Done');
	return
end

% Parameter file
DATA = LoadParameters([path '/' basename '.xml']);
disp(['... loaded parameter file ''' basename '.xml''']);

% Event file(s)
DATA.events.time = [];
DATA.events.description = {};
eventFiles = dir([path '/' basename '.*.evt']);
if ~isempty(eventFiles),
	for i = 1:length(eventFiles),
		events = LoadEvents([path '/' eventFiles(i).name]);
		if isempty(events.time), continue; end
		DATA.events.time = [DATA.events.time ; events.time];
		DATA.events.description = {DATA.events.description{:} events.description{:}}';
		disp(['... loaded event file ''' eventFiles(i).name '''']);
	end
else
	disp('... (no event file found)');
end

% Position file
DATA.positions = [];
if exist([path '/' basename '.pos']),
	DATA.positions = LoadPositions([path '/' basename '.pos'],DATA.rates.video);
	disp(['... loaded positions file ''' basename '.pos''']);
elseif exist([path '/' basename '.whl']),
	DATA.positions = LoadPositions([path '/' basename '.whl'],DATA.rates.video);
	disp(['... loaded positions file ''' basename '.whl''']);
elseif exist([path '/' basename '.whl']),
	DATA.positions = LoadPositions([path '/' basename '.mqa'],DATA.rates.video);
	disp(['... loaded positions file ''' basename '.mqa''']);
else
	disp('... (no positions file found)');
end

% Spike files
DATA.spikes = [];
for i = 1:DATA.spikeGroups.nGroups,
	filename = [path '/' basename '.' int2str(i) '.clu'];
	if exist(filename,'file'),
		try
			DATA.spikes = [DATA.spikes;LoadSpikeTimes(filename,DATA.rates.wideband)];
			disp(['... loaded spike files ''' basename '.' int2str(i) '.clu''']);
		catch
			disp(['... (could not load spike files ''' basename '.' int2str(i) '.clu'')']);
		end
	else
		filename = [path '/' basename '.clu.' int2str(i)];
		if exist(filename,'file'),
			try
				DATA.spikes = [DATA.spikes;LoadSpikeTimes(filename,DATA.rates.wideband)];
				disp(['... loaded spike files ''' basename '.clu.' int2str(i) '''']);
			catch
				disp(['... (could not load spike files ''' basename '.clu.' int2str(i) ''')']);
			end
		end
	end
end
if isempty(DATA.spikes),
	disp('... (no spike files found)');
end

% This is updated only once the files have been properly loaded
DATA.session.basename = basename;
DATA.session.path = path;

disp('Done');

