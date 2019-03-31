% lin_downsample        downsample with a MA
% 
% y = lin_downsample( x, dsf )
% 
% x     matrix (columns downsampled)
% dsf   downsampling factor (integer); negative number indices number of
%           samples in the output; vector indicates 

function y = lin_downsample( x, dsf )

y = [];
if ~isnumeric( x )
    return
end
if isvector( x )
    x = x( : );
end
n = size( x, 1 );
if length( dsf ) > 1
    didx = dsf;
else
    if dsf > 1 % downsample
        didx = round( dsf / 2 : dsf : n );
        MA = dsf;
    elseif dsf < -1
        ns = abs( dsf );
        didx = round( [ fliplr( n / 2 : -n/ns : 1 ) n / 2 + n / ns : n/ns : n ] );
        MA = ceil( mean( diff( didx ) ) );
    end
end
if ~exist( 'didx', 'var' )
    return
end
didx( 1 ) = max( didx( 1 ), 1 );
didx( end ) = min( didx( end ), n );
win = ones( MA, 1 ) / MA;
xf = firfilt( x, win );
y = xf( didx, : );

return