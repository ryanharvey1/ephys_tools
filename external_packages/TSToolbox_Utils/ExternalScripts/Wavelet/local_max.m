% LOCAL_MAX     detect all local extrema in a matrix
%
% [ idx, vals ] = local_max( x, mode )
%
% x         data matrix
% mode      (optional) string. 'ext', 'min'; default - maxima only
%
% idx       indices of local maxima
% vals      values
%
% note: if X is a matrix, this routine returns inds as a 2-column matrix
% (index of extremum, column index) and vals is a matched vector

% 15-jan-13 ES

function [ idx, vals ] = local_max( x, mode )

if nargin<2 || isempty(mode), mode = ''; end
x = colvec( x );
[ m n ] = size( x );
d2 = diff( sign( diff( x ) ) );
switch mode
case 'ext'
    [ row col ] = find( abs( d2 ) > 1 );
case 'min'
    [ row col ] = find( d2 > 1 );
otherwise
    [ row col ] = find( d2 < -1 );
end
row = row + 1;
if n == 1
    idx = row;
else
    idx = [ row col ];
end
vals = x( row + ( col - 1 ) * m );

return