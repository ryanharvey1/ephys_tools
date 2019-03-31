function [X,Y] = LoadPosition(pfn, varargin)

% [X,Y] = LoadPosition(pfn, parameters)
%
% INPUTS:
%	pfn is an ascii filename of position data of the form 
%		timestamp x y
%  pfn may also be a 'pvd' file, which has columns
%     usec_timestamp x y velocity direction
%     note that timestamps for this type of file are in microseconds
%     as opposed to the more standard 100 microsecond timestamp units
%     a file with the 'pvd' extension will have its timestamps divided by 
%     100 and the last two columns will be ignored.
% OUTPUTS:
%	X, Y tsArray data structures
%
% PARAMETERS:
%	ContinuousData 1/0 -- (default=0) if full allows use of ctsd structure

% ADR 1998 
% version L4.0
% status: PROMOTED
% Altered 4/2004 by David Euston so that it can load .pvd as well as .pascii
% files

%--------------------
% PARAMETERS

ContinuousData = 0;
Extract_varargin;

%-------------------------------------------
% Load file into t, x, and y

[fp, errmsg] = fopen(pfn, 'rt');
if (fp == -1)
  error(errmsg)
end

ReadHeader(fp);

[p fn ext] = fileparts(pfn);

if strcmpi(ext, '.pvd')
   disp('Reading PVD file format')
   A = fscanf(fp, '%f', [5,inf]); 
	t = A(1,:)' / 100;  % divide by 100 to convert from microseconds back to cheetah timestamps
	x = A(2,:)';
	y = A(3,:)';
else
	A = fscanf(fp, '%f', [3,inf]); 
	t = A(1,:)';
	x = A(2,:)';
	y = A(3,:)';
end;
fclose(fp);

%--------------------
% Write to standard format
%--------------------

if ContinuousData
   dt = mean(diff(t));
   X = ctsd(t(1), dt, x);
   Y = ctsd(t(1), dt, y);

   disp('WARNING (LoadPosition):');
   disp([ '        mean dt = ', num2str(mean(diff(t)))]);
   disp([ '        stddev = ', num2str(std(diff(t)))]);
else
   X = tsd(t, x);
   Y = tsd(t, y);
end