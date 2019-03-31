function WriteT(fn, ts)

% WriteT  Writes out an NSMA t-file
%
% WriteT(fn, ts)
% 
% INPUTS:
%   fn = filename to write to
%   ts = ts object to write out
% OUTPUTS:
%   (none)
%
% ADR 1998, version U1.0, last modified '98 by ADR

% status: PROMOTED


[fp,msg] = fopen(fn, 'wb', 'b');
if fp == -1; error(msg); end

WriteHeader(fp, 'T-file written by matlab');

fwrite(fp, Data(ts), 'uint32');

fclose(fp);