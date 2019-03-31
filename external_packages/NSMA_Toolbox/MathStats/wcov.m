function sigma = wcov(x, w)

% wcov  Returns weighted covariance matrix computed for pairs of columns of a matrix
%
% sigma = wcov(x, w)
%
% INPUTS:
%       x - nS x nD matrix
%       w - weights (must be nS x 1)
% OUTPUTS:
%       sigma - weighted covariance matrix (nD x nD) by col's of x (just like cov(x))
%
% ADR 1999, version L4.0, last modified '99 by ADR

% status: UNTESTED


%StatusWarning('UNTESTED', 'wcov');

if size(w) ~= size(x)
   error('Size mismatch in wmean');
end

[n,m] = size(x);

xc = (x - repmat(wmean(x, w), n, 1));
for iD = 1:m
   for jD = 1:m
      sigma(iD,jD) = (xc(:,iD)' * (w .* xc(:,jD))) / sum(w);
   end
end
