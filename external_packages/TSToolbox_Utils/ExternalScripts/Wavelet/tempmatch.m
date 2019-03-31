% TEMPMATCH         template matching for time series
%
% [ c, t, p ] = tempmatch( w, x, th )
%
% w         template (must be a vector)
% x         time series (can be matrix, works on columns)
% th        {0.9i} threshold: p-value (real) or corr (imag)
%
% c         correlation coefficient
% t         decision (1/0: match/mismatch)
% p         p-value
%
% for each point in data x, compute the correlation with template w
% then carry out an optional t-test (takes time..) to define matches

% 27-feb-12 ES

% revisions
% 06-jan-13 considerable speed-up using fft

function [ c, t, p ] = tempmatch( w, x, th )

c = [];
t = [];
p = [];

%------------------------------------------------------------%
% initialize
%------------------------------------------------------------%
nargs = nargin;
if nargs < 2 || isempty( w ) || isempty( x )
    error( 'missing arguments' )
end
if nargs < 3 || isempty( th ), th = 0.9j; end
[ n m ] = size( x );
w = w( : );
nw = length( w );
if ~isa( x, 'double' )
    x = double( x );
end
if n == 1
    x = x';
    n = m;
    m = 1;
end

%------------------------------------------------------------%
% compute the correlation with the template
%------------------------------------------------------------%
% prepare the filters
w = w - mean( w );
sw = sqrt( sum( w .* w ) / length( w ) );
w = w / sw;
awin = ones( nw, 1 ) / nw;
% determine nfft
nb = max( length( awin ), length( w ) );
fftflops = [ 18 59 138 303 660 1441 3150 6875 14952 32373 69762 ...
    149647 319644 680105 1441974 3047619 6422736 13500637 28311786 59244791];
nf = 2 .^ ( 1 : 20 );
validset = find( nf > ( nb - 1 ) );
nf = nf( validset );
fftflops = fftflops( validset );
L = nf - ( nb - 1 );
[ dum, ind ] = min( ceil( n ./ L) .* fftflops );
nfft = nf( ind );
L = L( ind );
% fft the filters
B1 = fft( awin, nfft );
B2 = fft( w( nw : -1 : 1 ), nfft );
B1 = B1( :, ones( 1, m ) );
B2 = B2( :, ones( 1, m ) );
% fft the signal blockwise
y1 = zeros( [ n m ] ); % awin, x
y2 = zeros( [ n m ] ); % awin, x^2
y3 = zeros( [ n m ] ); % w( inv ), x
x2 = x .* x;
istart = 1;
while istart <= n
    iend = min( istart + L - 1, n );
    yend = min( n, istart + nfft - 1 );
    if ( iend - istart ) == 0
        X = x( istart( ones( nfft, 1 ) ), : );
        X2 = x2( istart( ones( nfft, 1 ) ), : );
    else
        X = fft( x( istart : iend, : ), nfft );
        X2 = fft( x2( istart : iend, : ), nfft );
    end
    Y1 = ifft( X .* B1 );
    y1( istart : yend, : ) = y1( istart : yend, : ) + Y1( 1 : ( yend - istart + 1 ), : );
    Y2 = ifft( X2 .* B1 );
    y2( istart : yend, : ) = y2( istart : yend, : ) + Y2( 1 : ( yend - istart + 1 ), : );
    Y3 = ifft( X .* B2 );
    y3( istart : yend, : ) = y3( istart : yend, : ) + Y3( 1 : ( yend - istart + 1 ), : );
    istart = istart + L;
end
% compute the correlation
sx = sqrt( y2 - y1 .^ 2 );
c = y3 ./ sx / nw;
c = [ c( nw : n, : ); zeros( nw - 1, m ) ];
c( isinf( c ) ) = 0;

%------------------------------------------------------------%
% compute the p-values (from the correlation coefficients )
%------------------------------------------------------------%
if nargout >= 2
    dof = nw - 3;
    if nw > 50
        tstat = atanh( c ) .* sqrt( dof );
    else
        tstat = c .* sqrt( ( dof ) ./ ( 1 - c.^2 ) );
    end
    p = tcdf( -abs( tstat ), ( dof ) ); % two sided
    p( isinf( tstat ) ) = 0;
    if isreal( th )
        t = p < th;
    else
        t = c > imag( th );
    end
end

return

% EOF

% example:
x = rand( 1000, 2 ); 
w = x( 521 : 600, 2 ) + rand( 80, 1 );
[ c, t, p ] = tempmatch( w, x, 1e-5 );
find( t )
