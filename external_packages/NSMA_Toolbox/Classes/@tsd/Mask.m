function mTSD = Mask(tsd, varargin)

% tsd/Mask  Sets values in tsd object which are NOT within specified time periods to NaN
%
% mTSD = Mask(tsd, TrialPairs....)
%
% INPUTS:
%    tsd = tsd object
%    TrialPairs = pairs of start/end times (can be matrices of n x 2)
% OUTPUTS:
%    mTSD = tsd object with times *not* in TrialPairs set to NaN
%
% ADR 1998, version L4.0, last modified '98 by ADR

% RELEASED as part of MClust 2.0
% See standard disclaimer in ../Contents.m

% Unwrap trial pairs
MaskOFF = [StartTime(tsd)-1];
MaskON = []; 
for iTP = 1:length(varargin)
   curMask = varargin{iTP};
   MaskOFF = cat(1, MaskOFF, curMask(:,2));
   MaskON = cat(1,MaskON, curMask(:,1));
end
MaskON = cat(1, MaskON, EndTime(tsd)+1);
nTransitions = length(MaskON);
MaskON = sort(MaskON);
MaskOFF = sort(MaskOFF);

% Construction output tsd
mTSD = tsd;

% Now implement mask
for iT = 1:nTransitions
   f = find(mTSD.t > MaskOFF(iT) & mTSD.t < MaskON(iT));
   mTSD.data(f) = NaN;
end