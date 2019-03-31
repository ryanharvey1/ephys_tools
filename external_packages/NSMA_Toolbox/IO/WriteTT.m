function WriteTT(fn, TT)

% WriteTT  Writes an NSMA TT file
%
% WriteTT(fn, TT)
%
% INPUTS:
%   fn = name of .tt file to write
%   TT = tsd structure where data = nSpikes x nSamplesPerSpike x nTrodes
% OUTPUTS:
%   (none)
%
% uses mex file WriteTT0 to do the main write
%
% ADR 1998, version L1.0, last modified '98 by ADR

% status UNDER CONSTRUCTION


%StatusWarning('UNDER CONSTRUCTION', 'WriteTT');

t = Range(TT, 'ts');
wv = Data(TT);
if length(t) > 0
	WriteTT0(fn, t, wv);
else
   warning(['Number of Spikes = 0! No file ' fn ' was written!']);
end