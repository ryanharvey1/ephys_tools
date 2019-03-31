function R = overlap_AB(A,B)

% overlap_AB  Computes matrix of overlaps between pairs of columns in matrices A and B
%
% R = overlap_AB(A,B)
%
% INPUTS:
%       A,B - rectangular matrices with SAME NUMBER OF ROWS (= samples)
% OUTPUTS:
%       R = matrix of overlap coefficients - R_ij = overlap between each pair of columns A(:,i) and B(:,j)
%
% If  A = nR x nC1 and B = nR x nC2,  R will be nC1 x nC2 matrix of overlap coefficents in the range [0..1]
%
% PL 2001, last modified 4/23/01 by PL


[nRa, nCa] = size(A);
[nRb, nCb] = size(B);
if nRa ~= nRb
    error(' Matrices A and B must have same number of ROWS (samples)!' );
end

nrmdA = A./repmat(sqrt(sum(A.*A)),nRa,1); 
nrmdB = B./repmat(sqrt(sum(B.*B)),nRb,1);

R = (nrmdA' * nrmdB);
