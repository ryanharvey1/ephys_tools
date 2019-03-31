function sw = SpikeWidth(WV)
%
% sw = SpikeWidth(WV)
%
% INPUTS:
%     WV = tsd of nSpikes x nTrodes (4) x nChannels
%
% OUTPUTS:
%       sw = nSpikes x nTrodes (4) x spike width
%
% ALGO:
%       just time of valley - time of peak 
%       note: if valley preceeds peak, then sw is negative.

% ADR 1998
% version L5.1
% status PROMOTED

% v 5.0 -- completely different from 4.0 see new desc.
% v 5.1 -- vectorized code

WVD = Data(WV);

[peak, ipeak] = max(WVD, [], 3);
[vlly, ivlly] = min(WVD, [], 3);
sw = ivlly - ipeak;

