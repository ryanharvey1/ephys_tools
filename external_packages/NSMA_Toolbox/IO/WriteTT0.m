% WriteTT0  (MEX FILE) Writes an NSMA TT file
%
% WriteTT0(fn, t, wv)
%
% INPUTS: 
%   fn = name of .tt file to write
%   t = nSpikes x 1 timestamps
%   wv = nSpikes x nTrodes x nChannels array of waveforms
% OUTPUTS:
%   (none)
% 
% NOTE: Do not call this routine directly, use WriteTT (q.v.)
%
% ADR 1998, version L4.0, last modified '98 by ADR

% status: see C code