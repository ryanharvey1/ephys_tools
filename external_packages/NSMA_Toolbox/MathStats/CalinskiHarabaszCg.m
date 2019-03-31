function [Cg, B, W] = CalinskiHarabaszCg(mu, x, z)

% CalinskiHarabaszCg  Measure for quality of cluster separation
% 
% [Cg, B, W] = CalinskiHarabaszCg(mu, x, z)
%
% INPUTS:
%       mu - centers of clusters (nG=#clusters x nD=#dimensions) 
%       x - data (nS=#data points x nD=#dimensions)
%       z - assignments of data points to clusters (nS=#data points x 1)
% OUTPUTS:
%       Cg - C(g)=measure of separation of clusters
%       B - between groups dispersion
%       W - within groups dispersion
%
% ALGO
%     C(g) = (trace B / (nG - 1)) / (trace W / (nS - nG) )
%
% ADR 1999, version L4.0, last modified '99 by ADR

% status: PROMOTED


[nG, nD] = size(mu);
[nS, nD] = size(x);

% calculate W
W = zeros(nD, nD);
for iG = 1:nG
   z0 = (z == iG);
   for iS = 1:nS
      W = W + (z0(iS) .* (x(iS, :) - mu(iG, :)))' * (x(iS, :) - mu(iG, :));
   end
end

mx = mean(x);
% calculate B
B = zeros(nD, nD);
for iG = 1:nG
   B = B + sum(z == iG) * (mu(iG, :) - mx)' * (mu(iG, :) - mx);
end

Cg = (trace(B) / (nG - 1)) / (trace (W ) / (nS - nG));


   