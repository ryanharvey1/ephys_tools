% MY_SPECTRUM           Welch spectrum for multiple signals.
%
% call                  [ PXX, F ] = MY_SPECTRUM( X, NFFT, FS, WINDOW, NOVERLAP, DFLAG )
%
% gets                  X           matrix; trials in columns
%                       NFFT        {128}; length of segments to be averaged
%                       FS          {1000}; sampling rate in Hz
%                       WINDOW      vector for a window; otherwise hanning
%                                     of length WINDOW or {NFFT}
%                       NOVERLAP    samples overlap between segments in single column {NFFT/2}
%                       DFLAG       detrending mode: {'linear'}, 'none', or 'mean'
%
% returns               PXX         column-by-column auto spectrum
%                       F           frequency vector (matches PXX dimensions)
%
% does                  applies the Welch method of spectral estimation to
%                       each column (signal, trial) of X: divide signal
%                       into equal (overlapping) segments. window, fft, and
%                       power spectrum each segment. then average the
%                       spectra.
%
% calls                 nothing.

% 18-nov-12 ES

function [ Pxx, f, Psd ] = my_spectrum( X, nfft, Fs, window, noverlap, dflag )

% SS,SUM->SD
stdsumsq = inline( 'sqrt( max( (n*x2 - abs(x).^2)/(n-1), zeros(size(x)) ) )'...
    , 'x2', 'x', 'n' );

% arguments
nargs = nargin;
if nargs < 1 || isempty( X ), error( '1 argument' ), end
[ n m ] = size( X );                                            % data points x trials
if nargs < 2 || isempty( nfft ), nfft = min( n, 128 ); end
if numel( nfft ) > 1
    error( 'scalar nfft required' )
end
[ g e ] = log2( abs( nfft ) );
if g ~= 0.5
    nfft = 2 ^ e;
end
if nargs < 3 || isempty( Fs ) || numel( Fs ) ~= 1
    Fs = 1000;                                                  % 1 Khz
end
if nargs < 4 || isempty( window )
    window = hanning( nfft ); 
end
if numel( window ) == 1
    window = hanning( window ); 
end
window = window( : );
nwind = length( window );
window = window / sum( window );
if nargs < 5 || isempty( noverlap ) || ( noverlap > nfft / 2 )
    noverlap = ceil( length( window ) / 2 );                    % 50% overlap
end
if nargs < 6 || isempty( dflag )
    dflag = 'linear';
end

% prepare
if n < nwind                                                    % zero-pad x if length less than the window length
    X( n + 1 : nwind, : ) = 0;  
    n = nwind;
end
k = fix( ( n - noverlap ) / ( nwind - noverlap ) );             % number of windows
Pxx = zeros( nfft, m );                                         % initialize
Pss = Pxx;
index = 1 : nwind;

% accumulate sum
for i = 1 : k
    if strcmp( dflag, 'linear' )
        xd = detrend( X( index, : ) );
    elseif strcmp( dflag, 'none' )
        xd = X( index, : );
    else
        xd = detrend( X( index, : ), 0 );
    end
    xw = ( window * ones( 1, m ) ) .* xd;                       % apply smoothing window
    xx = fft( xw, nfft );                                       % fft columns
    Pxxi = abs( xx ) .^ 2;                                      % power
    Pxx = Pxx + Pxxi;                                           % accumulate segments
    Pss = Pss + Pxxi.^2;                                        % Sum of squares
    index = index + ( nwind - noverlap );                       % advance index
end

% output
if ~any( any( imag( X ) ~= 0 ) ),                               % if X not complex
    if rem( nfft, 2 )                                           % nfft odd
        select = [ 1 : ( nfft + 1 ) / 2 ];
    else
        select = [ 1 : nfft / 2 + 1 ];
    end
else
    select = 1 : nfft;
end

% compute SD from SS
if nargout > 2
    Psd = stdsumsq( Pss, Pxx, k );
    Psd = Psd( select, : ) / ( k * norm( window ) ^ 2 );            % select and normalize
end
Pxx = Pxx( select, : ) / ( k * norm( window ) ^ 2 );            % select and normalize
if nargout > 1
    f = ( select - 1 )' * Fs / nfft;                            % frequency vector
end

return

% EOF