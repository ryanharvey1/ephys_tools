% plotSpectrogram       plot an array of spectrograms on a single fig
% 
% fig = plotSpectrogram( x, y, z, USF, plotmode, zscaling, yscaling, nCycles )
%
% x, y, z       time/phase; frequency; power (or any other measure)
% USF           {1}; up-sampling factor
% plotmode      {'imagesc'} or 'contourf'
% zscaling      {'linear'}, 'log', or 'db'
% yscaling      {'linear'}, 'log'
% nCycles       {1.5} for phasograms
%
% can be called with empty x and/or y, but not with empty z
%
% calls         imupsample, alines, remrnd

% 13-mar-13 ES

% revisions
% 24-mar-13 moved the text to the bottom left panel
% to do
% (1) add the marginals 
% (2) support NaN

function [ fig ah zout ] = plotSpectrogram( x, y, z, USF, plotmode, zscaling, yscaling, nCycles, verbose )

sepColor = [ 1 1 1 ] * 1; % white
if nargout >= 3
    out = 1;
else
    out = 0;
    zout = [];
end

nargs = nargin;
if nargs < 3 || isempty( z )
    error( 'missing/empty argument' )
end
if isempty( x )
    x = 1 : size( z, 1 );
end
if isempty( y )
    y = 1 : size( z, 2 );
end
if length( x ) == size( z, 2 ) && length( y ) == size( z, 1 )
    z = permute( z, [ 2 1 3 ] );
end
if nargs < 4 || isempty( USF )
    if nargs > 4 && ~isempty( plotmode ) && isequal( plotmode, 'imagesc' )
        USF = 1; % default for image
    else
        USF = 10; % default for contourf
    end
end
if nargs < 5 || isempty( plotmode )
    %plotmode = 'imagesc';
    plotmode = 'contourf';
end
switch lower( plotmode )
    case {'image', 'imagesc' }
        plotmode = 'imagesc';
    case { 'contour', 'contourf' }
        plotmode = 'contourf';
end
if nargs < 6 || isempty( zscaling )
    zscaling = 'linear';
end
if nargs < 7 || isempty( yscaling )
    yscaling = 'linear';
end
if nargs < 8 || isempty( nCycles )
    nCycles = 1.5; % phase cycles (for phasograms)
end
if nargs < 9 || isempty( verbose )
    verbose = 0;
end

x = x( : ).';
if all( x >= -pi & x <= pi ) || all( x >= 0 & x <= 2 * pi )
    istime = 0;
    xlab = 'Phase [rad]';
else
    istime = 1;
    xlab = 'Time';
end
if ~istime
    x = [ x - 2 * pi  x x + 2 * pi ];
end

if strcmp( yscaling, 'log' )
    ytick = 2.^[ ceil( log2( y( 1 ) ) ) : floor( log2( y( end ) ) ) ];
    y = log10( y );
end
ylab = 'Frequency [Hz]';

n = size( z, 3 );
if n == 1
    newplot
    ah = gca;
    fig = gcf;
    klab = 1;
else
    nrows = ceil( n / ceil( sqrt( n ) ) );
    ncols = ceil( sqrt( n ) );
    [ ah fig ] = tilefig( nrows, ncols, 1, 0.85, 'bottomright' );
    klab = ( nrows - 1 ) * ncols + 1;
end

verb( sprintf( '%s: plotting # ', upper( mfilename ) ), -verbose )
for k = 1 : n
    verb( sprintf( '%d ', k ), -verbose )
    subplot( ah( k ) )
    if istime
        zk = z( :, :, k )';
    else
        zk = repmat( z( :, :, k )', [ 1 3 ] );
    end
    zk = local_scale( zk, zscaling );
    switch plotmode
        case 'imagesc'
            [ x1 y1 z1 ] = imupsample( x, y, zk, USF );
            imagesc( x1, y1, z1 );
            axis xy,
        case 'contourf'
            [ ignore ch ] = contourf( x, y, zk, USF ^ 2 );
            set( ch, 'linestyle', 'none' );
    end
    if out
        zout( :, :, k ) = zk;
    end
    title( sprintf( '%d: %0.3g', k, max( zk( : ) ) ) )
    if ~istime
        xlim( [ -1 1 ] * nCycles * pi )
        alines( -pi : pi : pi, 'x', 'color', sepColor, 'linestyle', '--', 'linewidth', 0.5 );
        set( ah( k ), 'xtick', ( -2 * pi : pi : 2 * pi )...
            , 'xticklabel', remrnd( ( -2 * pi : pi : 2 * pi ), 0.01 ) );
    end
    set( ah( k ), 'box', 'off', 'tickdir', 'out' )
    if k == klab
        ylabel( ylab )
        if n <= 6
            xlabel( xlab )
        end
        if strcmp( yscaling, 'log' )
            set( ah( k ), 'ytick', log10( ytick ), 'yticklabel',  ytick )
        end
        %set( ah( k ), 'xticklabel', '' )
    else
        set( ah( k ), 'xticklabel', '', 'yticklabel', '' )
    end
end
for k = ( n + 1 ) : size( ah, 1 )
    subplot( ah( k ) )
    axis off
end
colormap( myjet )
verb( sprintf( 'done!' ), verbose )

return

function z = local_scale( mat, scaling )
%bias = eps;
%bias = 0;
bias = NaN; % illegal hack for preventing visualizing 0's as NaNs:
if ismember( scaling, { 'log', 'db' } ) && isnan( bias )
    if sum( mat( : ) > 0 )
        mat( mat <= 0 ) = min( mat( mat > 0 ) );
    end
    bias = 0;
end
switch scaling
    case 'log'
        z = log10( abs( mat ) + bias );
    case 'db'
        z = 20 * log10( abs( mat ) + bias );
    case 'linear'
        z = mat;        
end
return

% EOF
