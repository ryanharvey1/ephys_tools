% MY_XCORR          compute the column-wise cross-correlation between real matrices.
%
% call              CC = MY_XCORR( X, Y, MLAG, NORMF )
%
% gets              X           matrix, trials in columns
%                   Y           {X}; another matrix of same dimensions; if
%                                   left empty, computes auto-correlation
%                   MLAG        {size(X,2)}; maximal lag
%                   NORMF       {0}, if 1, scale by autocorrelation at tau 0; 
%                                   if -1, compute cross-covariances   
%
% does              column-wise cross-correlation.
%
% NOTE:             if Y lags after X, peak is in negative lags
%                   (i.e. Y is the REFERENCE signal)
%
% calls             nothing.

% 23-feb-05 ES

% revisions
% 12-mar-05 default MLAG changed to full length
% 31-mar-05 IA normf option added 
% 01-apr-11 default MLAG changed to size( X, 2 ) - 1 
% 26-feb-12 negative normf computes in cross-covariance

function [ cc, lags ] = my_xcorr( x, y, MLAG, normf )

% arguments
nargs = nargin;
if nargs < 1 || isempty( x ), error( 'need input!' ), end
if nargs < 2 || isempty( y ), y = x; end
if nargs < 3, MLAG = []; end
if nargs < 4, normf = 0; end
if ~isequal( size( x ), size( y ) ), error( 'X, Y must be equal size' ); end
if MLAG < 0 || MLAG ~= round( MLAG ), error( 'MLAG must be non-negative scalar' ); end

% preparation
n = size( x, 1 );
if isempty( MLAG ), MLAG = n; end
MLAG = min( n - 1, MLAG );
lags = -MLAG : MLAG;
if normf < 0
    x = x - ones( n, 1 ) * mean( x );
    y = y - ones( n, 1 ) * mean( y );
end

% computation
nfft = 2 ^ nextpow2( 2 * n - 1 );
xx = fft( x, nfft );
yy = fft( y, nfft );
cc = real( ifft( xx .* conj( yy ) ) );
cc = [ cc( nfft - MLAG + 1 : nfft, : ); cc( 1 : MLAG + 1, : ) ];

if normf
    % Normalizes the sequence so that the auto-correlations
    % at zero lag are identically 1.0.
    if ~isequal(x,y),
        % Compute autocorrelations at zero lag
        cxx0 = sum(abs(x).^2);
        cyy0 = sum(abs(y).^2);
        scale = sqrt(cxx0.*cyy0); % geometric mean (as in xcorr)
        cc = cc./repmat(scale,size(cc,1),1);
    else
        % Autocorrelation case, simply normalize by c[0]
        cc = cc./repmat(cc(MLAG+1,:),size(cc,1),1);
    end
end

return

% example: equivalent to xcorr for a single trial  
x = rand( 100, 1 );
y = [ rand( 10, 1 ); x( 1 : 90, 1 ) ];
figure, plot( xcorr( x, y, 100 ), 'r' ), hold on, plot( my_xcorr( x, y ), ':b' )
% peak is in tau_vec(91) = -10