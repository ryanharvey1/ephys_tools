function m = wmean(x, w)

% wmean  Returns weighted means for each column of a matrix
%
% m = wmean(x, w)
%
% INPUTS:
%       x - matrix (nS x nD)
%       w - weights (must be nS x 1)
% OUTPUTS:
%       m - 1 x nD vector of weighted means by col (just like mean(x))
%
% ADR 1999, version L4.0, last modified '99 by ADR 

% status: UNTESTED


% StatusWarning('UNTESTED', 'wmean');

if size(w) ~= size(x)
   error('Size mismatch in wmean');
end

[nS, nD] = size(x);

for iD = 1:nD
   m(iD) = sum(w .* x(:,iD)) / sum(w);
end
