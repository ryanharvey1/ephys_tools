function SE = ReadSE_tsd(fn)

% ReadSE  Reads an NSMA SE file
%
% SE = ReadSE(fn)
%
% INPUTS:
%   fn = .SE file
% OUTPUTS:
%   SE = tsd Structure where data = nSpikes x nSamplesPerSpike x nTrodes
%
% Uses mex file ReadSE (PL '99) to do the main read
%
% ADR 2000, version L1.0, last modified '00 by ADR

% RELEASED as part of MClust 2.0
% See standard disclaimer in Contents.m


[t, wv] = ReadSE(fn);
SE = tsd(t, wv);