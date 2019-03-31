function L = SpkTrainValuation(S,f,r)

% L = SpkTrainValuation(S,f,t,r)
% 
% computes log-likelihood of spike train S with intensity function f.
% INPUTS:
%     S: spike timings during test epoch
%     f: a tsd describing the predicted intensity function during test epoch
%     r: average firing rate during training epoch
% 
% OUTPUT:
%     L: log-likelihood
%
% Dependencies: TStoolbox

% Adrien Peyrache, 2014 (following Harris, 2004)

t   = Range(f);
f   = Data(f);
dt  = median(diff(t));
ft  = tsd(t,f/r);
ft  = Restrict(ft,S);
d   = Data(ft);
intF    = sum(f-r).*(dt); %express it in sec!
d       = sum(log(d(d>0)));

L = d-intF;