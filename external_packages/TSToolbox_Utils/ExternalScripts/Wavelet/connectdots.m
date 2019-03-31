% connectdots           linear interpolation and extrapolation
%
% [ xc tc ] = connectdots( t, x )
%
% as name says, plus
% NaNs at the edges are extrapolated 

% 20-aug-13 ES

% revisions
% 21-aug-13 edges extrapolated

% NOTE:
% NaNs within the data are ignored (i.e. NaNs are output)
% should modify to be skipped


function [ xc tc ] = connectdots( t, x, g )

nargs = nargin;
if nargs < 2 || isempty( x ) || isempty( t )
    return
end
if nargs < 3 || isempty( g )
    g = 0;
end

n = length( t );
if length( x ) ~= n
    return
end
if ~issorted( t )
    [ t sidx ] = sort( t );
    t = round( t );
    x = x( sidx );
end

tc = ( t( 1 ) : t( end ) )';
m = length( tc );
nx = isnan( x );
xc = zeros( m, 1  );
for i = 1 : n - 1
    if i == 1 && nx( i ) && n > 2
        dy = x( i + 2 ) - x( i + 1 );
        dx = t( i + 2 ) - t( i + 1 );
        ni = t( i + 1 ) + 1 - t( i );
        x0 = x( i + 1 ) - diff( dy / dx * [ 1 ni ] );
        idx = ( t( i ) : t( i + 1 ) ) - t( 1 ) + 1;
    elseif i == ( n - 1 ) && nx( i + 1 ) && n > 1
        dy = x( i ) - x( i - 1 );
        dx = t( i ) - t( i - 1 );
        ni = t( i + 1 ) - t( i );
        x0 = x( i );
        idx = ( t( i ) : ( t( i + 1 ) - 1 ) ) - t( 1 ) + 1;
    else
        dy = x( i + 1 ) - x( i );
        dx = t( i + 1 ) - t( i );
        x0 = x( i );
        ni = dx;
        idx = ( t( i ) : ( t( i + 1 ) - 1 ) ) - t( 1 ) + 1;
    end
    xm = ones( ni, 1 ) * x0 + ( 0 : ni - 1 )' * dy / dx;
    xc( idx ) = xm;
end
if nx( n )
    xc( m ) = 2 * xc( m - 1 ) - xc( m - 2 );
else
    xc( m ) = x( n );
end

if logical( g )
    newplot
    line( t, x, 'marker', '.', 'color', [ 0 0 1 ] );
    line( tc, xc, 'linestyle', '--', 'color', [ 1 0 0 ] )
end

return

% EOF