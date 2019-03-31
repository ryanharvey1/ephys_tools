% SINDETECT     automatic detection of sine waves in an otherwise random signal
%
% CALL          [ pidx, aidx, trains, meandurs, si, zi, clu ] = sindetect( x, minAmp, minCC, minDC, minDur, rmvMed )
% 
% GETS          x           signal
%               minAmp      threshold for noise {1}
%               minCC       threshold for sine/zap match {[0.95 0.95]}
%               minDC       minimal instantaenous duty cycle; {0.5} (lower
%                               value -> more segmentation)
%               minDur      minimal segment duration {2} samples
%               rmvMed      {1}; by default, median is removed; flag this
%                               off for blockwise processing
%
% RETURNS       pidx		the peak times
%               aidx		the allocation of peaks to trains
%               trainidx	start and end times of each train
%               meandurs	the mean peak-to-peak intervals per train
%               scc         CC of piece-wise (segment) fit to sine
%               zcc         CC of global (train) fit to linear zap
%               clu         classification results: 1, sines; 2, zap; 0, neither
% 
% CALLS         parseSchmitt, calc_pearson, rankcols, alines

% ALGORITHM
% -detect segments above threshold (minAmp)
% -disard short ones (minDur)
% -merge adjacent segments into 'trains' (minDC)
% -check for each 'train' whether it is composed of (non/identical) sine waves (minCC)
% -classify into pure or variable-duration sines (1/2)

% 25-jan-13 ES

% revisions
% 27-jan-13 use rank correlation for robust support of various chirps
%           remove edges during regression
% 28-jan-13 (1) use a schimtt trigger (dual TH) detection mechanism
%           (2) compute duty-cycle from min (not mean) flanking segments
% 04-feb-13 minCycles implemented (false positive rejection)

function [ pidx, aidx, trains, meandurs, si, zi, clu ] = sindetect( x, minAmp, minCC, minDC, minDur, minCycles, rmvMed, graphics )

% initialize
pidx = [];
aidx = [];
trains = [];
meandurs = [];
si = [];
zi = [];
clu = [];

% arguments
nargs = nargin;
if nargs < 2 || isempty( minAmp )
    minAmp = 1;
end
if length( minAmp ) == 1
    minAmp = [ minAmp minAmp / 2 ];
end
if nargs < 3 || isempty( minCC )
    minCC = 0.95;
end
if length( minCC ) < 2
    minCC = minCC( 1 ) * [ 1 1 ];
end
if nargs < 4 || isempty( minDC )
    minDC = 0.5;
end
if nargs < 5 || isempty( minDur )
    minDur = 2;
end
if nargs < 6 || isempty( minCycles )
    minCycles = [ 2 10 ];
end
if length( minCycles ) == 1
    minCycles = minCycles * [ 1 1 ];
end
if nargs < 7 || isempty( rmvMed )
    rmvMed = 0;
end
if nargs < 8 || isempty( graphics )
    graphics = 0;
end

% prepare
x = x( : );
if rmvMed
    x = x - median( x );
end

% detect TH crossing segments
seg = parseSchmitt( x, minAmp( 1 ), minAmp( 2 ) );

% remove very short segments
durs = diff( seg, [], 2 ) + 1;
rmv = durs < minDur;
seg( rmv, : ) = [];
durs( rmv, : ) = [];

% partition into sequences, remove non-participating segments
m = size( seg, 1 );
if m < 2
    return
end
aidx = zeros( m, 1 );
gaps = seg( 2 : m, 1 ) - seg( 1 : m - 1, 2 );
mdurs = min( [ durs( 1 : m - 1 ) durs( 2 : m ) ], [], 2 );
dc = gaps ./ ( gaps + mdurs );
borders = [ 0; find( dc > minDC ); size( seg, 1 ) ];
b = length( borders );
trains = zeros( b - 1, 2 );
j = 0;
for i = 1 : ( b - 1 )
    if ( borders( i + 1 ) - borders( i ) ) == 1
        continue
    end
    j = j + 1;
    trains( j, : ) = [ seg( borders( i ) + 1, 1 ) seg( borders( i + 1 ), 2 ) ];
    aidx( borders( i ) + 1 : borders( i + 1 ) ) = j;
end
trains( j + 1 : end, : ) = [];
rmv = aidx == 0;
seg( rmv, : ) = [];
durs( rmv ) = [];
aidx( rmv ) = [];
m = size( seg, 1 );
if ~m
    return
end

% compute statistics for each remaining segment
pidx = zeros( m, 1 );
for i = 1 : m
	[ val idx ] = max( x( seg( i, 1 ) : seg( i, 2 ) ) );
	pidx( i, : ) = seg( i, 1 ) + idx - 1;
end

% check each train (sequence of segments) for fit with sine waves
ntrains = size( trains, 1 );
meandurs = zeros( ntrains, 1 );
si = zeros( ntrains, 1 );
zi = zeros( ntrains, 1 );
nsegs = zeros( ntrains, 1 );
for ti = 1 : ntrains
    idx = pidx( aidx == ti );
    if idx( 1 ) == trains( ti, 1 )
        idx( 1 ) = [];
    end
    if idx( end ) == trains( ti, 2 )
        idx( end ) = [];
    end
    nsegs( ti ) = length( idx ) - 1;
    if nsegs( ti ) <= 0
        continue
    end
    durs = diff( idx );
    meandurs( ti ) = mean( durs );
    durs( end ) = durs( end ) + 1;
    fs = 1 ./ durs;   
    x0 = x( idx( 1 ) : idx( end ) );
    nx = length( x0 );
    sinX = zeros( nx, 1 ); 
    cdurs = [ 0; cumsum( durs ) ];
    for i = 1 : length( durs ), 
        t = ( 0 : 1 : durs( i ) - 1 )';
        sinX( ( cdurs( i ) + 1 ) : cdurs( i + 1 ) ) = cos( 2 * pi * fs( i ) * t );
    end
    si( ti ) = calc_pearson( x0, sinX );
    if nsegs( ti ) > 2
        zi( ti ) = calc_pearson( ( 1 : length( durs ) )', rankcols( fs ) );
    end
    % plot( x0 ), hold on, plot( ( sinX + 1 ) / 2 * max( x0 ), 'r' )
end

% classify
clu = zeros( ntrains, 1 );
clu( si >= minCC( 1 ) & nsegs >= minCycles( 1 ) ) = 1; % any set of sine waves
clu( si >= minCC( 1 ) & abs( zi ) >= minCC( 2 ) & nsegs >= minCycles( 2 ) ) = 2; % any monotonous chirp

% plot
if graphics
    newplot
    plot( x )
    alines( trains( :, 1 ), 'x', 'color', [ 0 0.7 0 ] ); 
    alines( trains( :, 2 ), 'x', 'color', [ 1 0 0 ] );
    alines( minAmp, 'y', 'color', [ 1 1 1 ] * 0.5 );
    title( sprintf( '%d/%d/%d unclassified/sine/chirp'...
        , sum( clu == 0 ) ...
        , sum( clu == 1 ) ...
        , sum( clu == 2 ) ...
        ) );
end

return

% EOF
