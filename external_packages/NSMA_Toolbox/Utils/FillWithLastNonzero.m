function Y = FillWithLastNonzero(X)

% FillWithLastNonzero  Replaces all zero elements in a vector with the closest previous non-zero value
%
% Y = FillWithLastNonzero(X)
%
% INPUTS:
%       X - 1D array
% OUTPUTS:
%       Y - same as X, but all zero elements changed to closest previous non-zero value
%
% Note: currently uses for-loop.  There's got to be a better way.
%
% ADR 1998, version v4.0, last modified '98 by ADR

% Status: PROMOTED


Y = X;
lastX = NaN;
for iX = 1:length(X)
   if X(iX) 
      lastX = X(iX);
   else
      Y(iX) = lastX;
   end
end
