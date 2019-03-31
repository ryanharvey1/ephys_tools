% firFreqs          frequencies from fir filter

% 04-apr-13 ES

function freqs = firFreqs( hBP, Fs )

if isa( hBP, 'numeric' )
    hBP{ 1 } = hBP;
end
nfBins = length( hBP );
T = zeros( nfBins, 1 );
for i = 1 : nfBins
    ah = abs( hBP{ i } );
    minval = max( ah ) / 2;
    iStart = max( find( ah > minval, 1, 'first' ) - 1, 1 );
    iEnd = min( find( ah > minval, 1, 'last' ) + 1, length( ah ) );
    h = hBP{ i }( iStart : iEnd );
    v = diff( local_max( h, 'min' ) );
    if isempty( v )
        hh = hBP{ i };
        [ ign hT ] = min( hh( length( hh ) / 2 : end ) );
        T( i ) = 2 * hT;
    else
        T( i ) = mean( v );
    end
end
freqs = floor( Fs./T );

return

% EOF

