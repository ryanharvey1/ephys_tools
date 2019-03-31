% fft_upsample       integer upsample a signal in the frequency domain.
%
% y = fft_upsample( x, f, edges )
%
% f         {2} upsampling ratio
% edges     {0} flag
%
% does          
% works on matrix columns
% piecewise fft interpolation WITHOUT explicit filtering.
% default upsampling ratio is 2
% edges are taken care of by flipping
%
% calls         nothing
%
% notes         
% 1. X = Y( 1 : F : end, : ), within roundoff errors.
% 2. NaNs are supported (expanded to NaNs)
%
% see also lin_expand

% 08-jul-11 ES

function [ y xi ] = fft_upsample( x, f, edges )

nargs = nargin;
if nargs < 1, error( '1 argument' ), end
if nargs < 2 || isempty( f ), f = 2; end
if nargs < 3 || isempty( edges ), edges = 0; end

f = round( f );
m0 = size( x, 1 );
if edges
    x = [ flipud( x ); x; flipud( x ) ];
end
[ m n ] = size( x );
if m == floor( m / 2 ) * 2
    odd = 0;
else
    odd = 1;
    x = [ x; zeros( 1, n ) ];
    m = m + 1;
end
nans = isnan( x );
nans = nans( 1 : m, : );
x( nans ) = 0;
mf = m * f;
xfft = fft( x );
z = zeros( mf, n );
idx = [ 1 : m/2 mf - m/2 + 1 : mf ];
z( idx, : ) = xfft( 1 : m, : ) * f;
y = real( ifft( z ) );
nans = reshape( repmat( reshape( nans, [ 1 m n ] ), f, 1 ), [ m * f n ] );
y( nans ) = NaN;
if edges
    rmv = [ 1 : round( m / 3 * f ) round( 2 * m / 3 * f + 1 ) : m * f ];
    y( rmv, : ) = [];
end
if odd
    rmv = ( m0 * f + 1 ) : size( y, 1 );
    y( rmv, : ) = [];
end

if nargout > 1
    xi = ( 1 : 1/f : ( size( x, 1 ) + 1 - 1/f ) )';
end

return

% expanding NaNs by indexing:
v = logical( [ 0 0 0 1 1; 0 0 0 0 1 ]' );
v = reshape( v, [ 1 5 2 ] ); vv=repmat( v, 3, 1 ); reshape( vv, [ 15 2 ] ) 
% here, m = 5 rows, n = 2 columns, ratio f = 3

x = rand( 1228, 1 );
for f = 1 : 20; cc( f ) = calc_pearson( fft_upsample( x, f ), lin_expand( x, f ) ); end
subplot( 3, 1, 1 ), plot( x ), xlim( [ 1 length( x ) ] ), 
subplot( 3, 1, 2 ), plot( fft_upsample( x, 10 ) ), xlim( [ 1 10 * ( length( x ) - 1 ) ] ), 
subplot( 3, 1, 3 ), plot( lin_expand( x, 10 ) ), xlim( [ 1 10 * ( length( x ) - 1 ) ] )
