% MA_RMS            moving average RMS of a signal.
%
%                   Y = MA_RMS( X, M, RMSFLAG )
%
%                   RMS computed in a moving average window of M samples:
%
%                   RMS(n) = sqrt( 1 / M * sum( ( x(n-k) )^2 ) )
%
%                   where the summation is over k=-M/2:M/2 (odd M)
%
%                   1. if input is matrix, works on columns
%                   2. NaNs are taken care of
%                   3. zero phase lag
%
% calls             firfilt

% 19-nov-12 ES

function y = ma_rms( x, m, rmsflag )

nargs = nargin;
if nargs < 2 || isempty( m ), m = size( x, 1 ); end % time independent
if nargs < 3 || isempty( rmsflag ), rmsflag = 1; end % 0: MS

if m > size( x, 1 )
    y = x; 
    return
end

nans = isnan( x );
if sum( nans( : ) )
    flag = 1;
    x( nans ) = 0;
else
    flag = 0;
end
x2 = x .^ 2;
win = ones( m, 1 ) / m;
y = firfilt( x2, win );
if rmsflag
    y = sqrt( y );
end
if flag
    y( nans ) = NaN;
end

return