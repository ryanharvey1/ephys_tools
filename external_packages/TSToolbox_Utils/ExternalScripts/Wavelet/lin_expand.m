% lin_expand        linear expansion without any filtering
%
% y = lin_expand( x, f )
%
% see also: fft_upsample, interp1, linint

% 03-jan-13 ES

function y = lin_expand( x, f )

f = round( f( 1 ) );
sx = size( x );
sy = sx;
sy( 1 ) = sy( 1 ) * f;
x = x( : );
y = repmat( x, [ 1 f ] )';
y = reshape( y, sy );

return

% EOF