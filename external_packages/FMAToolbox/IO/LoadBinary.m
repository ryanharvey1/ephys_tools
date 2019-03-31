%LoadBinary - Load data from a multiplexed binary file.
%
%  USAGE
%
%    data = LoadBinary(filename,<options>)
%
%    filename       file to read
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'duration'    duration to read (in s) (default = Inf)
%     'frequency'   sampling rate (in Hz) (default = 20kHz)
%     'start'       position to start reading (in s) (default = 0)
%     'nChannels'   number of data channels in the file (default = 1)
%     'channels'    channels to read (default = all)
%     'precision'   sample precision (default = 'int16')
%     'skip'        number of bytes to skip after each value is read
%                   (default = 0)
%    =========================================================================

% Copyright (C) 2004-2010 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function data = LoadBinary(filename,varargin)

% Default values
start = 0;
nChannels = 1;
precision = 'int16';
skip = 0;
duration = Inf;
frequency = 20000;
channels = [];

if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+3) ' is not a property (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'duration',
			duration = varargin{i+1};
			if ~isdscalar(duration,'>=0'),
				error('Incorrect value for property ''duration'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		case 'frequency',
			frequency = varargin{i+1};
			if ~isdscalar(frequency,'>0'),
				error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		case 'start',
			start = varargin{i+1};
			if ~isdscalar(start),
				error('Incorrect value for property ''start'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
			if start < 0, start = 0; end
		case 'nchannels',
			nChannels = varargin{i+1};
			if ~isiscalar(nChannels,'>0'),
				error('Incorrect value for property ''nChannels'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		case 'channels',
			channels = varargin{i+1};
			if ~isivector(channels,'>=0'),
				error('Incorrect value for property ''channels'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		case 'precision',
			precision = varargin{i+1};
			if ~isstring(precision),
				error('Incorrect value for property ''precision'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		case 'skip',
			skip = varargin{i+1};
			if ~isiscalar(skip,'>=0'),
				error('Incorrect value for property ''skip'' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help LoadBinary">LoadBinary</a>'' for details).']);
	end
end

% By default, load all channels
if isempty(channels),
	channels = 1:nChannels;
end

% Check consistency between channel IDs and number of channels
if any(channels>nChannels),
	error('Cannot load specified channels (listed channel IDs inconsistent with total number of channels).');
end

sizeInBytes = 0;
switch precision,
	case {'uchar','unsigned char','schar','signed char','int8','integer*1','uint8','integer*1'},
		sizeInBytes = 1;
	case {'int16','integer*2','uint16','integer*2'},
		sizeInBytes = 2;
	case {'int32','integer*4','uint32','integer*4','single','real*4','float32','real*4'},
		sizeInBytes = 4;
	case {'int64','integer*8','uint64','integer*8','double','real*8','float64','real*8'},
		sizeInBytes = 8;
end

if ~exist(filename),
	error(['File ''' filename ''' not found.']);
end
f = fopen(filename,'r');
if f == -1,
	error(['Cannot read ' filename ' (insufficient access rights?).']);
end


% Position file index for reading
start = floor(start*frequency)*nChannels*sizeInBytes;
status = fseek(f,start,'bof');
if status ~= 0,
	fclose(f);
	error('Could not start reading (possible reasons include trying to read past the end of the file).');
end

if isinf(duration),
	% Determine number of samples when duration is 'inf'
	fileStart = ftell(f);
	status = fseek(f,0,'eof');
	if status ~= 0,
		fclose(f);
		error('Error reading the data file (possible reasons include trying to read past the end of the file).');
	end
	fileStop = ftell(f);
	nSamplesPerChannel = (fileStop-fileStart)/nChannels/sizeInBytes;
	duration = nSamplesPerChannel/frequency;
	frewind(f);
	status = fseek(f,start,'bof');
	if status ~= 0,
		fclose(f);
		error('Could not start reading (possible reasons include trying to read past the end of the file).');
	end
else
	nSamplesPerChannel = round(frequency*duration);
	if nSamplesPerChannel ~= frequency*duration,
%  		disp(['Warning: rounding duration (' num2str(duration,'%.15g') ' -> ' num2str(nSamplesPerChannel/frequency,'%.15g') ')']);
		duration = nSamplesPerChannel/frequency;
	end
end

% For large amounts of data, read chunk by chunk

maxSamplesPerChunk = 100000;
nSamples = nChannels*nSamplesPerChannel;
if nSamples > maxSamplesPerChunk,
	% Determine chunk duration and number of chunks
	nSamplesPerChunk = floor(maxSamplesPerChunk/nChannels)*nChannels;
	durationPerChunk = nSamplesPerChunk/frequency/nChannels;
	nChunks = floor(duration/durationPerChunk);
	% Preallocate memory
	data = zeros(nSamplesPerChannel,length(channels));
	% Read all chunks
	i = 1;
	for j = 1:nChunks,
		d = LoadBinaryChunk(f,'frequency',frequency,'nChannels',nChannels,'channels',channels,'duration',durationPerChunk,'skip',skip);
		[m,n] = size(d);
		if m == 0, break; end
		data(i:i+m-1,:) = d;
		i = i+m;
	end
	% If the data size is not a multiple of the chunk size, read the remainder
	remainder = duration - nChunks*durationPerChunk;
	if remainder ~= 0,
		d = LoadBinaryChunk(f,'frequency',frequency,'nChannels',nChannels,'channels',channels,'duration',remainder,'skip',skip);
		[m,n] = size(d);
		if m ~= 0,
			data(i:i+m-1,:) = d;
		end
	end
else
	if skip ~= 0,
		data = fread(f,[nChannels frequency*duration],precision,skip);
	else
		data = fread(f,[nChannels frequency*duration],precision);
	end
	data=data';

	if isempty(data),
		warning('No data read (trying to read past file end?)');
	elseif ~isempty(channels),
		data = data(:,channels);
	end
end

fclose(f);
