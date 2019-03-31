function L = SpkTrainLogLikelihood(q,f)

% L = SpkTrainValuation(S,f,r)
% 
% computes log-likelihood of spike train S with intensity function f.
% INPUTS:
%     q: binned spike train
%     f: predicted intensity function (in number of spikes)
% 
% OUTPUT:
%     L: log-likelihood

% Adrien Peyrache, 2014 (following Harris, 2004)

L = -f+q.*log(f);
L = sum(L(f>0));
if isinf(L) || isnan(L)
   L =  0;
end