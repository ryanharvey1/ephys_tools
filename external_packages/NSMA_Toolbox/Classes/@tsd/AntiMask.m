function mTSD = AntiMask(tsd, varargin)

% tsd/AntiMask  Sets values in tsd object which are WITHIN specified time periods to NaN
%
% mTSD = AntiMask(tsd, TrialPairs....)
%
% INPUTS:
%    tsd = tsd object
%    TrialPairs = pairs of start/end times (can be matrices of n x 2)
% OUTPUTS:
%    mTSD = tsd object with times in TrialPairs set to NaN
%
% ADR 1998, version L4.0, last modified '98 by ADR

% RELEASED as part of MClust 2.0
% See standard disclaimer in ../Contents.m

MaskON = [];
MaskOFF = [];
% Unwrap trial pairs
for iTP = 1:length(varargin)
   curMask = varargin{iTP};
   MaskON = cat(1,MaskON, curMask(:,1));
   MaskOFF = cat(1, MaskOFF, curMask(:,2));
end
nTransitions = length(MaskON);

% Construction output tsd
mTSD = tsd;

% Now implement mask
for iT = 1:nTransitions
   f = find(mTSD.t > MaskON(iT) & mTSD.t < MaskOFF(iT));
   mTSD.data(f) = NaN;
end