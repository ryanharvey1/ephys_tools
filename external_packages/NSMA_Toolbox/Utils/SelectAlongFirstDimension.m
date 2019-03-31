function R = SelectAlongFirstDimension(IN, f)

% SelectAlongFirstDimension  Restricts data of n-dimensional array to specified indices in the 1st dimension
%
% R = SelectAlongFirstDimension(IN,f)
%
% INPUTS:
%       IN = any matrix input
%       f = selection indices (such as returned by find)
% OUTPUTS:
%       R = same type matrix as IN 
%
% equivalent of IN(f) for 1D, IN(f,:) for 2D, IN(f,:,:) for 3D, etc...
%
% ADR 1998, version L4.0, last modified '98 by ADR

% status PROMOTED


sz = size(IN);
dim1 = sz(1);
dimRest = sz(2:length(sz));

tmpMatrix = reshape(IN, [dim1 prod(dimRest)]);
tmpMatrix = tmpMatrix(f,:);

R = reshape(tmpMatrix, [length(f) dimRest]);