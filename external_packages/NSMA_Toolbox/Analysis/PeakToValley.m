function p2v = PeakToValley(WV)

% p2v = PeakToValley(WV)
%
% INPUTS:
%     WV = tsd of nSpikes x nTrodes (4) x nChannels
%
% OUTPUTS:
%       p2v = nSpikes x nTrodes (4) x spike width
%
% ALGO:
%       ratio of peak height:valley depth

% ADR 1998
% version L1.0
% status PROMOTED

WVD = Data(WV);
peak = max(WVD, [], 3);
vlly = min(WVD, [], 3);
S = warning;
warning off
p2v = abs(peak)./abs(vlly);
warning(S);
