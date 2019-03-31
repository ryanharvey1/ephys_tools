% psine             pulsed sine waveform
%
% x = psine( n, f )
% 
% a psine with (n,0) is a perfect sine wave
% a psine with (n,1) is a perfect pulse
% a psine with (n,0.5) is a n/2 pulse flanked by an n/2 sine cut into two halfs
%
% all are scaled to the 0-1 range

% 03-feb-13 ES

function x = psine( n, f )

x = [];

nargs = nargin;
if nargs < 2 || isempty( n ) || isempty( f )
    return
end
n = round( n( 1 ) );
f = f( 1 );
if n < 0 || ~isfinite( n ) || ~isfinite( f )
    return
end
if f < 0
    f = 0;
end
if f > 1
    f = 1;
end

n1 = floor( round( n * f ) / 2 ) * 2;
r = sin( 2 * pi * ( 1 : n1 )' / n1 - pi / 2 );
x = [ r( 1 : n1 / 2 ); ones( n - n1, 1 ); r( n1 / 2 + 1 : n1 ) ];
x = ( x + 1 ) / 2;

return