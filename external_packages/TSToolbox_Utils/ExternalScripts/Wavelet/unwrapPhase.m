% unwrapPhase        into discrete cycles
%
% receives
% phs       vector, phases 0-2*pi
% t0        center
%
% uphs      unwrapped phases
%
% does
% -looks for the trough (phase pi) closest to t0. defines this as the
% trough of cycle0
% -goes forward through the data, define a peak as a 0-phase crossing
% preceeded by a pi-phase crossing (and vice versa for troughs)
% -do the same backwards
% -define cycles between peaks


function [ uphs talign ] = unwrapPhase( phs, t0, graphics )

uphs = [];
talign = [];

nargs = nargin;
if nargs < 1 || isempty( phs )
    fprintf( '%s: phs missing\n', upper( mfilename )  )
    return
end
phs = phs( : );
n = length( phs );
if nargs < 2 || isempty( t0 )
    t0 = 1;
end
uphs = phs;
talign = t0;
if t0 > n || t0 < 1 || length( t0 ) > 1 || t0 ~= round( t0 ) 
    fprintf( '%s: t0 mismatch\n', upper( mfilename )  )
    return
end
if nargs < 3 || isempty( graphics )
    graphics = 0;
end


% detect the cycle0 trough
zc = phs - pi;
trf = find( zc( 1 : n - 1 ) <= 0 & zc( 2 : n ) >= 0 );
if isempty( trf )
    return
end
[ minval minidx ] = min( abs( t0 - trf ) );
trf0 = trf( minidx );
% find first peak before it
uzi = unwrap( phs );
pk0 = floor( uzi( trf0 ) / ( 2 * pi ) ) * 2 * pi;
%cyc0 = find( uzi < pk0, 1, 'last' ) + 1;
uphs = uzi - pk0;
talign = t0 - trf0; % positive number - shift the trigger back



% to reconstruct the cycle number (on a sample-by-sample basis):
% cycnum = floor( uphs / ( 2 * pi ) )
% the potential issue with that is wriggling, but that can be handled
% easily, see EOF

if graphics
    
    
    % detect the cycle0 trough
    zc = phs - pi;
    trf = find( zc( 1 : n - 1 ) <= 0 & zc( 2 : n ) >= 0 );
    [ minval minidx ] = min( abs( t0 - trf ) );
    trf0 = trf( minidx );
    % get the other relevant events
    uzi = unwrap( phs );
    trf1 = 2 * pi * ceil( min( uzi ) / ( 2 * pi ) ) + pi;
    tphs = trf1 : 2 * pi : max( uzi );
    pk1 = 2 * pi * ceil( min( uzi ) / ( 2 * pi ) );
    pks = pk1 : 2 * pi : max( uzi );
    % combine the tagged (0: peak; 1: trough) events and sort
    events = [ [ zeros( length( pks ), 1 ) pks( : ) ] ; [ ones( length( tphs ), 1 ) tphs( : ) ] ];
    events = sortrows( events, 2 );
    nevents = size( events, 1 );
    events = [ events NaN * ones( nevents, 1 ) ];
    % go forward/backwards and determine event times
    ev0 = [ 1 uzi( trf0 ) ];
    [ ign ev0idx ] = min( abs( events( :, 2 ) - ev0( 2 ) ) );
    events( ev0idx, 3 ) = trf0;
    for eidx = ( ev0idx +  1 ) : nevents
        t = find( uzi >= events( eidx, 2 ), 1, 'first' );
        events( eidx, 3 ) = t;
    end
    for eidx = ( ev0idx -  1 ) : -1 : 1
        t = find( uzi < ( events( eidx, 2 ) ), 1, 'last' ) + 1;
        events( eidx, 3 ) = t;
    end
    % prune and check errors:
    nans = isnan( events( :, 3 ) );
    events( nans, : ) = [];
    if sum( abs( diff( events( :, 1 ) ) ) ~= 1 ) > 0
        error( 'algorithmic issues..' )
    end
    
    % bias the phases
    sidx = events( events( :, 1 ) == 0, 3 ); % peaks
    if sidx( 1 ) > 1
        sidx = [ 1; sidx ];
    end
    eidx = [ sidx( 2 : end ) - 1; n ];
    cyc0 = find( trf0 > sidx & trf0 < eidx );
    ncyc = length( sidx );
    bias = [ -cyc0 + 1 : 0  1 : ncyc - cyc0 ]';
    %[ sidx eidx bias ]
    % now enumerate
    b = zeros( n, 1 );
    for i = 1 : ncyc
        b( sidx( i ) : eidx( i ) ) = bias( i );
    end
    % add to the phase
    uphs = phs + 2 * pi * b;
    % keep the offset
    talign = t0 - trf0; % positive number - shift the trigger back
    % e.g if ctrig=2, the trigger has to be moved 2 samples (Fout) back
    
    
    figure
    subplot( 2, 2, 1 )
    plot( phs )
    alines( 0 : pi : 2 * pi, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );
    alines( t0, 'x' );
    
    subplot( 2, 2, 2 )
    plot( uzi ); alines( pks, 'y', 'color', [ 1 0 0 ] ); alines( tphs, 'y', 'color', [ 0 0 1 ] );
    alines( events( events( :, 1 ) == 0, 3 ), 'x', 'color', [ 1 0 0 ] );
    alines( events( events( :, 1 ) == 1, 3 ), 'x', 'color', [ 0 0 1 ] );
    
end


return

% to reconstruct the cycle number (on a sample-by-sample basis):
% cycnum = floor( uphs / ( 2 * pi ) )
cycnum = floor( uz / ( 2 * pi ) );
% however, this may be non-monotonic due to wriggling. to make monotonic:
monotonic( cycnum )

% not perfect because accounts for wriggles only in the forward direction;
% so the full thing is:
cycnum = floor( uz / ( 2 * pi ) );
cychat = [ -flipud( monotonic( flipud( -cycnum( cycnum < 0 ) ) ) ); cycnum( cycnum == 0 ); monotonic( cycnum( cycnum > 0 ) ) ];
plot( cycnum, 'b' );
hold on, plot( cychat, 'r' );
alines( 0, 'y', 'linestyle', '--', 'color', [ 0 0 0 ] );

x0 = cycnum;
dx = [ 0; diff( x0 ) ];
neg = find( dx < 0 );
for i = 1 : length( neg )
    xidx = neg( i ) : length( x0 ); 
    x0( xidx ) = x0( xidx ) + 1;
end



% EOF

