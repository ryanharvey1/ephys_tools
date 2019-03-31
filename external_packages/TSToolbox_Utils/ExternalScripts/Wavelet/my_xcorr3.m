% MY_XCORR3         3 signal cross correlation by cross bispectrum.
%
% compute the cross correlation between 3 stationary signals
%
%   (1)     Rxyz(m,n) = E[ x(k)y(k+m)z(k+m+n) ]
%
% for the 3 signals case, the extension of WK theorem is
%
%   (2)     Sxyz( wk, wl ) = sum sum Rxyz( m,n ) * exp( -j*wk*m ) * exp( -j*wl*n )
%                              k = 0, 1, ... p-1
%                              l = 0, 1, ... p-1
%
% when the cross bispectrum is
%
%   (3)     Sxyz( wk, wl ) = X( wk ) * Y( wl ) * Z*(wk + wl ).
%
% NOTE
%
% for the 2 signal case, matlab defines cross correlation as 
%
%   (4)     Rxy(m) = E[ x(n+m)y(n ) ],
%
% whereas i use 
%
%   (5)     Rxy(m) = E[ x(n)y(n+m) ].
%
% this results in a lag reversed vector ( or a matrix flipped in both
% dimensions for the 3 signal case ).
%
% ALGORITHM
%
%   1. compute the cross bispectrum ( fft and multiply, (3) )
%   2. ifft the result ( WK, (2) )

% CONVERGENCE
% the direct method used here provides a consistent, approximately
% unbiased estimator of the cross bispectrum. methods to reduce variance
% include:
% 1. increasing the number of records (eg ifr xcorr3; see CALC_IFR_X)
% 2. smoothing in the frequency domain ( results in loss of resolution and
% bias )
% 3. increasing signal length (optimal if signal is stationary...)

% COMPLEXITY
% assume 3 signals of length n
% with desired offset m; define k = 2*m+1
% 1. the naive algorithm is O(n^3)
% 2. improved enumeration is O(n*k^2)
% 3. the fft aproach is:
%   3.1 indexing: O( n^2 )
%   3.2 fft:      3*n*log(n) = O( nlog(n) )
%   3.3 ifft:     2*n^2*log(n^2) = O( n^2log(n) )
%   therefore the fft approach is O( n^2log(n) ), which
%   is always better than (1) but not necessarily (2).
%
% if k^2 < nlog(n) => (2)
% else             => (3).
% 
% to summarize:
% for offsets that approach or surpass n/2, use fft;
% for small offsets, smart enumeration is faster.

% REFERENCES
% 1. matlab's XCORR
% 2. Nikias CL & Mendel JM, Signal processing with higher order spectra. IEEE
% signal processing magazine, july 1993.

% 09-nov-03 ES

% next
% 1. optimize by switching between algorithms as described above
% 2. still need to check fliplr of the result.

function [ cc, lags ] = my_xcorr3( x1, x2, x3, maxlag )

if nargin ~= 4
    error( '4 arguments' )
end

% get signal lengths
x1 = x1(:); 
x2 = x2(:); 
x3 = x3(:);
n = length( x1 );
n2 = length( x2 );
n3 = length( x3 );
if ~isequal( n, n2, n3 )
    error( 'signals must be of equal length' )
end
if ~isreal( [ x1; x2; x3 ] )
    error( 'signals must be real' )
end
% these requirements may be easily removed

% clip maxlag to data
nfft = 2^nextpow2( 2*n - 1 );
if maxlag > nfft
    maxlag = nfft - 1;
    warning( 'clipping lags to data length' )
end

% create indices for 3rd signal's fft
ci = 1 : nfft; 
ri = [ 0 : nfft - 1 ]'; 
sidx = ri( :, ones( 1, nfft ) ) + ci( ones( nfft, 1 ),: );
didx = [ 1 : nfft 1 : nfft - 1 ];
idx = didx( sidx );

% cross bispectrum
X1 = fft( x1, nfft );
X2 = fft( x2, nfft );
X3c = conj( fft( x3, nfft ) );
S123 = (X1 * X2.') .* X3c( idx );

% cc
ccomp = ifftn( S123 );

% force real
cc = real( ccomp );

% return only desired indices
lags = -maxlag : maxlag;
si = 1 : maxlag + 1;
ei = nfft - maxlag + 1 : nfft;
cc = [ cc( ei, ei ) cc( ei, si ); cc( si, ei ) cc( si, si ) ];

return
