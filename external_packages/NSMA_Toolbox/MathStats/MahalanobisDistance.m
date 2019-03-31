function M = MahalanobisDistance(mu, covar)

% MahalanobisDistance  Computes Mahalanobis distances between N-dimensional gaussian distributions
%
% M = MahalanobisDistance(mu, covar)
%
% INPUTS:
%       mu = array nG x nD of centers of gaussians
%       covar = array nD x nD x nG of covariance matrices
% OUTPUTS:
%       M = nG x nG array of Mahalanobis distances
%
% MahalanobisDistance is the part in the exponent exp(mah_dist)^2 of an N-dimensional gaussian.
% If x is a column vector (point in N-dim space), mu the mean column vector and Cinv 
%   the inverse of the covariance matrix, mah_dist = ((x-mu)'*Cinv*(x-mu))/2 .
%
% The definition of Mahalanobis distance between gaussians i and j =
%    ( (mu_j - mu_i)' * covar_i^-1 * (mu_j - mu_i))
%
% REFERENCE McLachlan and Basford (1988) Mixture Models: Inference and applications 
% to Clustering. Marcel Dekker Inc. New York and Basel.
%
% ADR 1999, version L4.0, last modified '99 by ADR

% status: PROMOTED


[nG,nD] = size(mu);
if [nD,nD,nG] ~= size(covar)
   error('MahalanobisDistance: Input size mismatch.');
end

M = zeros(nG, nG);

for iG = 1:nG
   sigma = covar(:,:,iG);   
   for jG = 1:nG
      M(iG, jG) = (mu(jG, :) - mu(iG, :)) * inv(sigma) * (mu(jG, :) - mu(iG, :))';
   end
end

M = sqrt(M);