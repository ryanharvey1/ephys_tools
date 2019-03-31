% MA_CORR           moving average correlation between multiple signals.
%
% CC = MA_CORR( X, M )
%
% pearson correlation computed in a moving average window of M samples
%
% 1. if input is matrix, works on columns
% 2. NaNs are NOT taken care of
% 3. zero phase lag
%
% calls             firfilt

% 09-sep-13 ES

% work on this...

function cc = ma_corr( x, m, vflag )
%m = 1000;

nargs = nargin;
if nargs < 1 || isempty( x ), return, end
[ rows cols ] = size( x );

if nargs < 2 || isempty( m ) || m > rows, m = rows; end % time independent
if nargs < 3 || isempty( vflag ), vflag = 0; end % 

% prepare filters, pairs, allocate memory..
win = ones( m, 1 );
awin = win / m;
pairs = make_ai( cols );
npairs = size( pairs, 1 );
cc = zeros( rows, npairs );
tim = zeros( 1, npairs );

% single column computations
if vflag
t0 = clock;
fprintf( '%d columns. preps...', cols )
end
x0 = x - firfilt( x, awin );
x2 = x0 .^ 2;
sd = sqrt( firfilt( x2, win ) );
if vflag
fprintf( '%0.3g sec. \tcomputing ', etime( clock, t0 ) )
end

% pair-wise computations
for i = 1 : npairs
    if vflag
    fprintf( '%d ', i )
    t0 = clock;
    end
    p1 = pairs( i, 1 );
    p2 = pairs( i, 2 );
    num = x0( :, p1 ) .* x0( :, p2 );
    num = firfilt( num, win );
    den = sd( :, p1 ) .* sd( :, p2 );
    cc( :, i ) = num ./ den;
    if vflag
    tim( i ) = etime( clock, t0 );
    end
    %fprintf( '%0.2g', tim( i ) );
end
if vflag
fprintf( 'done (%0.3g sec)\n', sum( tim ) )
end
%cc = firfilt( cc, awin );

return

% EOF

cc78 = zeros( ( size( x, 1 ) - m ), 1 );
for i = 1 : ( size( x, 1 ) - m )
    xx = x( i + [ 1 : m ], 7 : 8 );
    xm = bsxfun( @minus, xx, mean( xx ) );
    cc78( i ) = sum( xm( :, 1 ) .* xm( :,2 ) ) / sqrt( sum( xm( :, 1 ) .^ 2 ) * sum( xm( :, 2 ) .^ 2 ) );
    %calc_pearson( xx );
    % xm = bsxfun( @minus, xx, mean( xx ) );
    % sum( xm( :, 1 ) .* xm( :,2 ) ) / sqrt( sum( xm( :, 1 ) .^ 2 ) * sum( xm( :, 2 ) .^ 2 ) )
end
    
    
tic
cc0 = [];
for i = 1 : ( size( x, 1 ) - m )
    xx = x( i + [ 1 : m ], : );
    % xm = bsxfun( @minus, xx, mean( xx ) );
    % sum( xm( :, 1 ) .* xm( :,2 ) ) / sqrt( sum( xm( :, 1 ) .^ 2 ) * sum( xm( :, 2 ) .^ 2 ) )
    %cc0( i ) = calc_pearson( xx );
    mat = calc_pearson( xx ); 
    mat( diagidx( mat ) ) = 0; 
    cc0( i, : ) = squareform( mat );
end
toc

cc1 = cc( ( 1 : size( cc0, 1 ) ) + ceil( m / 2 ), : );

% computationally - this is an enormous squeeze. 
% e.g. for 8 channels, 1e5 samples, this takes just 0.7 sec
% in contrast, the brute-force approach takes about 200 sec

