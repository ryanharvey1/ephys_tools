% MYDETREND         remove DC and detrend columns of a matrix/3D array
%
% Y = mydetrend( Y, dflag )
%
% dflag     {1}; dc removal flag (if off, only detrends)

% 28-mar-11 ES

% revisions
% 28-apr-11 fast matrix implementation
% 03-sep-13 dflag added

function Y = mydetrend( Y, dflag )

[ M N ] = size( Y );
if nargin < 2 || isempty( dflag )
    dflag = 1;
end
if dflag
    dc = ones( 1, M );
    dc = dc / norm( dc );
    m = dc * Y;
    d = ( dc' * ones( 1, N ) ) .* ( ones( M, 1 ) * m );
    Y = Y - d;
end
tr = -( M - 1 ) / 2 : ( M - 1 ) / 2;
tr = tr / norm( tr );
m = tr * Y;
d = ( tr' * ones( 1, N ) ) .* ( ones( M, 1 ) * m );
Y = Y - d;

return

% EOF