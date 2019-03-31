% calcPhase         compute single-frequency phase
%
% [ phs, xf ] = calcPhase( x, fir, method, usf, graphics )
%
% ARGUMENTS:
% x             signal (vector or matrix). this can be:
%                   -the actual data (then any method can be applied)
%                   -a list of sample numbers corresponding to extrema (then 
%                       only 'peaks' and 'troughs' are valid). This mode is
%                       indicated by an imaginary input (complex(0,x)). if
%                       multiple segments are utilized then a second column
%                       indicates the segment indices (see local_max)
%               'troughs' can be applied). in the latter case, the 
% fir           {[]}; 3-element vector: [ fmin fmax Fs ] 
%                   Note that if usf > 1, then Fs is assumed to apply to
%                   the frequency of the upsampled signal.
%                   -Alternatively a vector, or a cell array of vectors
%                   (again, if usf > 1, the fir is applied to the upsampled
%                   data)
% method        {'hilbert'}, 'peaks', 'troughs', 'extrema', or 'wavelet'
% usf           {1}; upsampling factor; if integer >1, x is upsampled before 
%                   phase extraction 
% graphics      {0}; flag
%
% OUTPUT:
% phs           the phase, resolved 0-2pi. 0/2pi corresponds to peak, pi to trough
% xf            optional second output xf is the filtered version of x
%
% DOES
% -filters signal, then applies a hilbert transform to get the phases
% -if the method is 'peaks'/'troughs'/'extrema', applies a linear interpolation 
%   between the extrema instead of the hilbert transform. Thus phase is undetermined
%   before/after the first/last extermum
% -if the method is 'wavelet', then applies a wavelet transform (morlet
%   mother wavelet). in that case, fir must be a 3-element vector. 
%   In that case the phase is slightly distorted at the edges
%
% see also      eegPhase
%
% calls         fft_upsample, makefir, firfilt, local_max, wavelet, scale

% 10-feb-13 ES

% to do:
% -implement the 'extrema' method (assign 0 to peaks, pi to troughs, then
%       interpolate linearly between)

function [ phs, xf ] = calcPhase( x, fir, method, usf, graphics )

phs = [];
xf = [];

% arguments
nargs = nargin;
nout = nargout;
if nargs < 1 || isempty( x ), return; end
if nargs < 2 || isempty( fir ), fir = []; end
if nargs < 3 || isempty( method ), method = 'hilbert'; end
method = lower( method );
if nargs < 4 || isempty( usf ), usf = 1; end
usf = ceil( usf );
if nargs < 5 || isempty( graphics ), graphics = 0; end

% prepare/filter
if isvector( x )
    x = x( : );
end
if isequal( imag( x( :, 1 ) ), abs( x( :, 1 ) ) )
    issignal = 0;
else
    issignal = 1;
end
if issignal && usf ~= 1
    x = fft_upsample( x, usf );
end
[ nrows ncols ] = size( x );
if issignal
    if isempty( fir )
        return
    end
    switch method
        case 'wavelet'
            if length( fir ) ~= 3
                return
            end
            fMin = fir( 1 ) / sqrt( 2 );
            fMax = fir( 2 ) * sqrt( 2 );
            Fs = fir( 3 ) * usf;
            nbins = 2;
            dt = 1 / Fs;
            s0 = 1 / fMax;
            tMax = 1 / fMin;
            dj = log2( tMax/s0 ) / nbins;
            mother = 'MORLET';
        case { 'peaks', 'troughs', 'hilbert' }
            % apply filter/s
            if isa( fir, 'cell' )
                xf = x;
                for i = 1 : length( fir )
                    if ~isreal( fir{ i } )
                        if length( fir{ i } ) == 1
                            x0 = medfilt1( xf, imag( fir{ i } ) );
                        else
                            x0 = firfilt( xf, imag( fir{ i } ) );
                        end
                        xf = xf - x0;
                    else
                        xf = firfilt( xf, fir{ i } );
                    end
                end
            elseif isvector( fir )
                if length( fir ) == 3
                    fBP = fir( 1 : 2 );
                    Fs = fir( 3 );
                    fir = makefir( fBP, Fs * usf, [], 'bandpass' );
                elseif length( fir ) < 3
                    return
                end
                xf = firfilt( x, fir );
            end
        otherwise
            return
    end
elseif ~ismember( method, { 'peaks', 'troughs' } )
    return
end

% compute phase
switch method
    case { 'peaks', 'troughs' }
        if strcmp( method, 'peaks' )
            bias = 0;
            ext = 'max';
        else
            bias = pi;
            ext = 'min';
        end
        if issignal
            ts = local_max( xf, ext );
            cols = 1 : ncols;
        else
            if size( x, 2 ) == 2
                ts = abs( x );
            else
                tim = imag( x );
                tag = ones( size( tim, 1 ), 1 ) * [ 1 : ncols ];
                ts = [ tim( : ) tag( : ) ];
            end
            if usf > 1
                ts( :, 1 ) = ts( :, 1 ) * usf;
            end
            cols = unique( ts( :, 2 )' );
            ncols = length( cols );
            nrows = max( ts( :, 1 ) );
        end
        if size( ts, 2 ) == 1
            ts = [ ts ones( length( ts ), 1 ) ];
        end
        if issignal
            phs = NaN * ones( [ nrows ncols ] );
        else
            phs = [];
        end
        for col = cols
            t = ts( ts( :, 2 ) == col, 1 );
            n = length( t );
            ph = cumsum( 2 * pi * ones( n, 1 ) );
            ti = ( t( 1 ) : t( end ) )';
            phi = interp1( t, ph, ti, 'linear' );
            if issignal
                phs( ti, col ) = mod( phi, 2 * pi ) - bias;
            else
                phs = [ phs; phi col * ones( length( phi ), 1 ) ];
            end
        end
    case 'hilbert'
        phs = angle( hilbert( xf ) );
    case 'wavelet'
        if nout > 1 || graphics
            xf = zeros * ones( [ nrows ncols ] );
        end
        phs = zeros * ones( [ nrows ncols ] );
        for col = 1 : ncols
            xw = wavelet( x( :, col ), dt, 1, dj, s0, nbins, mother );
            if nout > 1 || graphics
                xf( :, col ) = real( xw( 2, : ) )';
            end
            phs( :, col ) = angle( xw( 2, : ) )';
        end
end
phs = mod( phs, 2 * pi );

% graphics
if graphics
    newplot
    i = 1;
    if issignal
        t = 1 : length( x( :, i ) );
        plot( t, scale( x( :, i ) ), 'b', t, scale( xf( :, i ) ), 'r', t, scale( phs( :, i ) ), 'k' )
        [ cc ll ] = my_xcorr( xf, x, length( x ), -1 ); 
        [ ign maxidx ] = max( cc ); 
        tstr = sprintf( '; max @ %d', ll( maxidx ) );
        legend( 'raw', 'filtered', 'phase' )
        alines( find( diff( sign( mod( phs, 2 * pi ) - pi ) ) ), 'x', 'color', [ 1 0 0 ], 'linestyle', '--' );
    else
        plot( phs( :, 1 ), 'k' );
        alines( find( [ 1; diff( phs( :, 2 ) ) ] ), 'x', 'color', [ 0 0 1 ] );
        ylabel( 'phase (rad)' )
        legend( 'phase', method )
        tstr = '';
    end
    xlabel( 'sample number' ) 
    title( sprintf( '%s%s', method, tstr ) )
end

return

% EOF

% to test and compare:
Fs = 20000;
t = ( 1 / Fs : 1 / Fs : 1/10 )';
x = sin( 2 * pi * 100 * t );
pos = [ 37         643        1884         461 ];

figure, phs1 = calcPhase( x, [ 50 150 Fs ], 'hilbert', [], 1 ); set( gcf, 'position', pos );
figure, phs2 = calcPhase( x, [ 50 150 Fs ], 'peaks', [], 1 ); set( gcf, 'position', pos );
figure, phs3 = calcPhase( x, [ 50 150 Fs ], 'troughs', [], 1 ); set( gcf, 'position', pos );
figure, phs4 = calcPhase( x, [ 50 150 Fs ], 'wavelet', [], 1 ); set( gcf, 'position', pos );
figure, phs5 = calcPhase( complex( 0, local_max( x, 'max' ) ), [ 50 150 Fs ], 'peaks', [], 1 ); set( gcf, 'position', pos );
figure, phs6 = calcPhase( complex( 0, local_max( x, 'min' ) ), [ 50 150 Fs ], 'troughs', [], 1 ); set( gcf, 'position', pos );
