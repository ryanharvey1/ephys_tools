function [PeakData, PeakNames,PeakPars] = feature_PEAK15to25(V, ttChannelValidity, Params)

% MClust
% [PeakData, PeakNames] = feature_PEAK15to25(V, ttChannelValidity)
% Calculate peak feature max value for each channel, using only samples
% 15-25.
%
% INPUTS
%    V = TT tsd
%    ttChannelValidity = nCh x 1 of booleans
%
% OUTPUTS
%    Data - nSpikes x nCh peak values
%    Names - "Peak15to25: Ch"
%

% ADR April 1998
% Modified by BMH Aug 2017 - use only inds from 15 to 25
% version M1.0
% See standard disclaimer in Contents.m



TTData = V.data();

[nSpikes, nCh, nSamp] = size(TTData);


f = find(ttChannelValidity);


PeakData = zeros(nSpikes, length(f));

PeakNames = cell(length(f), 1);
PeakPars = {};
PeakData = squeeze(max(TTData(:, f, 15:25), [], 3));

for iCh = 1:length(f)
   PeakNames{iCh} = ['Peak15to25: ' num2str(f(iCh))];
end
