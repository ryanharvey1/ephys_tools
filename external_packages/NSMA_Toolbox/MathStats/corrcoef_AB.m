function R = corrcoef_AB(A,B)

% corrcoef_AB  Computes matrix of correlations between pairs of columns in matrices A and B
%
% R = corrcoef_AB(A,B)
%
% INPUTS:
%       A,B - rectangular matrices with SAME NUMBER OF ROWS (= samples)
% OUTPUTS:
%       R = matrix of correlations - R_ij = correlation between each pair of columns A(:,i) and B(:,j)
%
% If  A = nR x nC1  and B = nR x nC2,  R will be nC1 x nC2 matrix of corrlation coefficents
% This function is a generalization of the matlab function corrcoef(A) which corresponds to corrcoef_AB(A,A)
%
% PL 2001, last modified 4/23/01 by PL


[nRa, nCa] = size(A);
[nRb, nCb] = size(B);
if nRa ~= nRb
    error(' Matrices A and B must have same number of ROWS (samples)!' );
end

stdA = (A-repmat(mean(A),nRa,1))./(repmat(std(A),nRa,1));
stdB = (B-repmat(mean(B),nRb,1))./(repmat(std(B),nRb,1));

R = (stdA' * stdB)/(nRa-1);
