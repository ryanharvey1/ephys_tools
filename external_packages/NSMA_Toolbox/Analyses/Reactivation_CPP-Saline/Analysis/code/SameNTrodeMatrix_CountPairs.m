function ncount = SameNTrodeMatrix_CountPairs(sameNTMatrix)
%  counts number of unique pairs in a same-NTrode-matrix, without double counting
%
% ncount = SameNTrodeMatrix_CountPairs(sameNTMatrix)
%
% returns number of pairs in a same-NTrode-matrix (which is 1 for a valid pair, zeros otherwise)
%
% PL Feb 2003

% count only upper triangle of square matrix:

[nr, nc] = size(sameNTMatrix);
if nr ~=nc
    error('sameNTMatrix must be a square matrix with zeros and ones');
end

ncount = 0;
for i=1:nr
   ncount = ncount + sum(diag(sameNTMatrix,i-1));
end