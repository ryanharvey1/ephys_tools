function y = logpoisspdf(x, lambda)

% logpoisspdf  Returns log(poisspdf(x,lambda)) [log(P(X=x)) where X~Poisson(lambda)]
%
% y = logpoiss(x, lambda)
%
% INPUTS:
%       x - value to find probability of in Poisson(lambda) distribution
%       lambda - parameter (mean) of Poisson distribution
% OUTPUTS:
%       y - log(P(X=x)) where X~Poisson(lambda)
%
% more efficient than poisspdf for very small probabilities
%
% ADR 1998, version v4.0, last modified '98 by ADR

% status: PROMOTED

if lambda == 0
   y = zeros(size(x));
   return;
end

y = zeros(size(x));
for iX = 1:length(x)
   y(iX) = -lambda + x(iX) * log(lambda) - sum(log(1:x(iX)));
end