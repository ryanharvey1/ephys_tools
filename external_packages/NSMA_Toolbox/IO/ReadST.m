function ST = ReadST(fn)

% ReadST  Reads an NSMA ST file
%
% ST = ReadST(fn)
%
% INPUTS:
%   fn = .ST file% OUTPUTS:%   ST = tsd structure where data = nSpikes x nSamplesPerSpike x nTrodes%% Uses mex file ReadST (PL '99) to do the main read%% ADR 2000, version L1.0, last modified '00 by ADR
% RELEASED as part of MClust 2.0% See standard disclaimer in Contents.m
[t, wv] = ReadST(fn);ST = tsd(t, wv);