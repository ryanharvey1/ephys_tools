% DFTFILT           digital filtering using FFT
%
% call              [ Y FR ] = DFTFILT( X, FS, FLP, FHP, BAND )
%
% gets              X       matrix, signals in columns
%                   Fs      sampling frequency
%                   FLP     lowpass cutoff
%                   FHP     highpass cutoff
%                   BAND    {1} for bandpass; 0 for bandstop
%
% returns           Y       matrix, filtered signals in columns
%                   FR      fraction of signal power in the filtered band
%
% calls             nothing
%
% algorithm         zeros all irrelevant FFT coefficients
%                   
% note              this is a very weak filter

% 07-apr-13 ES

function [ y fr ] = dftfilt( x, fs, flp, fhp, band )

nargs = nargin;
if nargs < 2 || isempty( x ) || isempty( fs ), error( 'missing input' ), end
if nargs < 5 || isempty( band ), band = 1; end
if any( size( x ) == 1 ), x = x( : ); end

N = size( x, 1 );
f = ( 0 : N - 1 ) / N * fs;

if nargs < 3 || ( isempty( flp ) && isempty( fhp ) )
    error( 'at least one band limit must be specified' )
elseif nargs == 3 || isempty( fhp )
    fzero = f > flp;                        % lowpass
elseif nargs >= 4 && isempty( flp ),
    fzero = f < fhp;                        % highpass
elseif nargs >= 4 && band == 1
    fzero = f < flp | f > fhp;              % bandpass
elseif nargs >= 4
    fzero = f > flp & f < fhp;              % bandstop
end

xfft = fft( x );
xfft( fzero, : ) = 0;
y = real( ifft( xfft ) );

if nargout > 1
    fr = ( var( x ) - var( x - y ) ) ./ var( x );
end

return

% EOF

