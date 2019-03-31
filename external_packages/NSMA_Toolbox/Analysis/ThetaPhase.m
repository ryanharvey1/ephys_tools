function [ph, thpeaks] = ThetaPhase(S, CRtheta, ThStart, ThEnd)
% ph = ThetaPhase(S, CRtheta, ThStart, ThEnd)
%
% Computes the theta phase of each spike for a group of cells 
%
% INPUTS:
% S:              cell array of ts objects containing the spike trains
% CRtheta:        EEG tsd object filtered for theta
% ThStart, Thend: arrays containig the start and stop times of the valid
%                 theta epochs
%
% OUTPUTS: 
% ph:             cell array of tsd abjects containing the  timestamp of each
%                 spike in Range(ph,'ts') and the phase S (in range [0,1] and the theta cycle number
%                 in a 2 column array in Data(ph)
% thpeaks:        The timestamps of each theta peak (the peak sample point of CRtheta tsd)
%
% batta 2000 (modified by PL Aug. 2000)
% status: alpha

t = Range(CRtheta, 'ts');
dd = Data(CRtheta);
dth = diff(dd);

dth1 = [0 dth'];
dth2 = [dth' 0];
clear dth;


peakindex = find(dth1 > 0 & dth2 < 0);
thpeaks = t(peakindex);
thpeaksAmplitudes = dd(peakindex);


clear t;
ph = cell(length(S),1);
for iC = 1:length(S)
%   SiC = Restrict(S{iC}, ThStart, ThEnd); 
%   s = Data(SiC);
  s = Data(S{iC});
  ph{iC} = zeros(size(s));
  pks = zeros(size(s));
  for j = 1:length(s)
    pk = binsearch_floor(thpeaks, s(j));
    if pk ~= length(thpeaks)
    ph{iC}(j) = (s(j) - thpeaks(pk)) / (thpeaks(pk+1) - thpeaks(pk));
    pks(j) = pk;
    end
  end
  ph{iC} = tsd( s, [ph{iC} pks]);
end


