function make_evt_file(event_ts,varargin)
% make_evt_file: for use with neuroscope

p = inputParser;
p.addParameter('eventtype','epochs'); %string for you that = event type you're saving
p.addParameter('numeventtypes',2); % you can have 3 diff types of events, start, peak, stop
p.addParameter('description',{'start','stop'}); 
p.addParameter('basePath',pwd);
p.parse(varargin{:});

eventtype = p.Results.eventtype;
numeventtypes = p.Results.numeventtypes;
description = p.Results.description;
basePath = p.Results.basePath;

parts = strsplit(basePath,filesep);
baseName = parts{end};

% Save as .evt for inspection
baseName = [baseName eventtype '.RO1.evt']; % you need the ROX bc neuroscope is buggy and uses this to parse files.

% Below is for ripples specifically
% Populate events.time field
lengthAll = size(event_ts,1)*numeventtypes;
events.time = zeros(1,lengthAll);

for i = 1:numeventtypes
    events.time(i:numeventtypes:lengthAll) = event_ts(i,:);
end

% Populate events.description field
events.description = cell(1,lengthAll);

for i = 1:numeventtypes
    events.description(i:numeventtypes:lengthAll) = description(i);
end

% Save .evt file for viewing in neuroscope - will save in your current directory
SaveEvents(fullfile(basePath,baseName),events) %Save and load into neuroscope along with corresponding LFP file
end

function SaveEvents(filename,events)

%SaveEvents - Write events to file.
%
%  USAGE
%
%    SaveEvents(filename,events)
%
%    filename            event file name
%    events              event data
%
%  SEE
%
%    See also NewEvents, LoadEvents, SaveRippleEvents.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if exist(filename), error('File already exists. Aborting.'); end

file = fopen(filename,'w');
if file == -1,
	error(['Cannot write to ' filename]);
end

for i = 1:length(events.time),
	fprintf(file,'%f\t%s\n',events.time(i)*1000,events.description{i}); % Convert to milliseconds
end

fclose(file);
end