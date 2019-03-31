function [ts_env,pow_env, PeaksIdx] = PowerEnvelope(ts,pow)
%        [ts_env,pow_env, PeaksIdx] = PowerEnvelope(ts,pow)
%
% calculate the envelope (= set of all local maxima) of the local power
% (= amplitude squared or smoothed version of amplitude squared) of a (filtered)
% eeg signal
%
% input: pow  ... local power (amplitude squared) N x 1 array
%        ts   ... corresponding timestamps        N x 1 array
%
% output: ts_env ... timestamps of envelope points (subset of ts) Ns x 1 array
%         pow_env ... envelope of local power pow  Ns x 1 array
%         PeaksIdx ... array of indizes of local maxima in pow
%
%  PL Dec. 99


[nr,nc]= size(pow);
if nc == 1
   % pow is a column vector; do nothing
elseif nr == 1
   pow = pow'; % pow is a row vector; transpose!
else
   error('input pow must be a column vector');
end

de = diff(pow);
de1 = [de; 0];   % backward derivative
de2 = [0; de];   % forward derivative
   
%finding peaks
PeaksIdx = find(de1 < 0 & de2 > 0);  % finds indices of local peaks
ts_env = ts(PeaksIdx);
pow_env = pow(PeaksIdx);
 