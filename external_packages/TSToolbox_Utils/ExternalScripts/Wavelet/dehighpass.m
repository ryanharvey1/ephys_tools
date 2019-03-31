% dehighpass        in a data-dependent manner
%
% xd = dehighpass( x, minAmp, minDur, f )
%
% default parameters:
% minAmp            0.01 of the range (max-median)
% minDur            3; [samples]
% f                 2; multiples of the instantaneous duration
% 
% does
% given a segment of data filtered with a high-pass filter (say 2-pole
% Butterworth), detects "pulses" (of thershold crossing) and compensates
% piecewise linearly for the filter's exponential decline
%
% parameters
% minAmp should be above the noise level (about 6mV for AmpliPex, so 0.02
% is fine). minDur can be 0, but if there is occasional spiky noise would
% be better to be 2 or 3 samples. f is non-trivial
%
% << how to decide on f >>
% -short f (say <1) would not group components of a train; since single-element
% trains are ignored, no correction would be made
% -long f may group multiple trains so the first element of a train may be
% too high
% -thus f should be as high as possible, but no higher than the minimum
% inter-train interval divided by the element duration
% -better to err on the high side
% -2 or 3 is usually a good compromise
% 
% see also          parseSchmitt (calls SchmittTrigger)
%                   parse (calls parsec mex file)
%                   connectdots
%                   firfilt

% 20-aug-13 ES

% revisions
% 21-aug-13 extapolation of edges

% NOTE
% edges are extrapolated LINEARLY, whereas the filtering is exponential
% it doesn't matter for the end, but the beginning of each segment is
% undercorrected
% ideally, an exponential model should be fitted to the dots (instead of
% connectdots.m); this is much more involved and may still be an
% underestimate since the initial part will be ignored anyhow

function xd = dehighpass( x, minAmp, minDur, f, g )

% parameters:
nargs = nargin;
if nargs < 1 || isempty( x )
    xd = [];
    return;
end
x = x( : );
xd = x;
if nargs < 2 || isempty( minAmp )
    minAmp = 0.01 * [ max( abs( x ) ) - median( x ) ];
end
if length( minAmp ) == 1
    minAmp = minAmp * [ 1 0.9 ];
end
if length( minAmp ) ~= 2 || minAmp( 2 ) >= minAmp( 1 )
    return
end
if nargs < 3 || isempty( minDur )
    minDur = 3; % [samples]
end
if nargs < 4 || isempty( f )
    f = 3;
end
if nargs < 5 || isempty( g )
    g = 0;
end

% partition into trains
mat = parseSchmitt( x, minAmp( 1 ), minAmp( 2 ) );      % events
T = diff( mat, 1, 2 ) + 1;                              % durations
rmv = T < minDur;
mat( rmv, : ) = []; 
T( rmv, : ) =   [];
if isempty( mat )
    return
end
m = size( mat, 1 );
t1 = mat( :, 1 ); 
t2 = mat( :, 2 );
dt = t1( 2 : m ) - t2( 1 : m - 1 );                     % gaps
mT = ( T( 1 : m - 1 ) + T( 2 : m ) ) / 2;
idx = find( dt < f * mT );
if isempty( idx )
    return
end
nat = parse( idx );
nat( :, 2 ) = nat( :, 2 ) + 1;
rmv = diff( nat, 1, 2 ) == 1;
nat( rmv, : ) = [];
n = size( nat, 1 );

for i = 1 : n

    % correct
    idx = nat( i, 1 ) : nat( i, 2 );
    k = length( idx );
    ttims = round( ( t2( idx( 1 : k - 1 ) ) + t1( idx( 2 : k ) ) ) / 2 );
    pad = ceil( [ T( idx( 1 ) ) T( idx( k ) ) ] / 2 );
    edges = [ t1( idx( 1 ) ) t2( idx( k ) ) ] + pad .* [ -1 1 ];
    edges( 1 ) = max( edges( 1 ), 1 );
    edges( 2 ) = min( edges( 2 ), length( x ) );
    tims = [ edges( 1 ); ttims; edges( 2 ) ];
    vals = [ NaN; x( ttims ); NaN ];
    lower = connectdots( tims, vals );
    w = round( mean( T( idx ) ) );
    win = ones( w, 1 ) / w;
    lf = firfilt( lower, win );
    xidx = tims( 1 ) : tims( end );
    xd( xidx ) = x( xidx ) - lf;

    % plot selected
    if g && i == g
        newplot
        subplot( 2, 1, 1 )
        line( xidx, x( xidx ), 'color', [ 0 0 0 ] )
        line( ttims, x( ttims ), 'marker', '.', 'color', [ 0 0 1 ] )
        line( xidx, lf, 'color', [ 0 1 0 ] )
        subplot( 2, 1, 2 )
        line( xidx, xd( xidx ) )
        line( xidx( [ 1 end ] ), [ 0 0 ], 'linestyle', '--', 'color', [ 0 0 0 ] )
    end

end

return

% EOF

figure, 
subplot( 2, 1, 1 ), plot( x ), hold on, plot( xd, 'r' )
subplot( 2, 1, 2 )
plot( x, xd, '.' ), sqaxis( gca, 1 ); xlabel( 'x' ), ylabel( 'xd' ), 
title( sprintf( 'CC=%0.2g; %0.3g%% modified', calc_pearson( x( : ), xd( : ) ), 1 - sum( x == xd ) / length( x ) ) )
