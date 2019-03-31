function [x,y,ci] = calcXYData(hObj,xlims)
% calcXYData -- Calculate X- and Y-Data.
%
%   [x,y,ci] = calcXYData(hObj,xlims) computes x- and y-data and (optional)
%   confidence intervals for the given the xlims (x-limits).
%
%   [x,y,ci] = calcXYData(hObj) looks to the data space of the axes to get the
%   x-limits.
%
%   If xlims is empty, then the empty is returned for x, y and ci.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2010/12/20 00:09:41 $

% If no x-limits are supplied, then try to get them from the data space attached
% to the axes.
if nargin < 2
    xlims = i_GetXLimFromAxesDataSpace( hObj );
end

% If we don't have any axes limits, then return empty for each argument.
if isempty( xlims )
    x  = zeros( 0, 1 );
    y  = zeros( 0, 1 );
    ci = zeros( 0, 2 );
    return
end

% Start the x-data but getting linearly spaced values between the x-limits.
x = linspace(xlims(2),xlims(1),hObj.Granularity).';

% Add real x values, if known, so interpolants will pass through data
ua = hObj.UserArgs;
if iscell(ua) && length(ua)>=1
    dshandle = ua{2};
    if ~isempty(dshandle) && ishandle(dshandle)
        if isa(hObj.fit,'cftool.fit')
            xdata = getxdata(dshandle,hObj.fit.outlier);
        else
            xdata = getxdata(dshandle,'');
        end
        x = sort([x(:); xdata(:)]);
    end
end

% If requested, try to compute confidence bounds
ci = [];
y = [];
ws = warning('off', 'all');
warningCleanup = onCleanup( @() warning( ws ) );
if isequal(hObj.ShowBounds,'on') && hObj.dfe>0
    try
        [ci,y] = predint(hObj.Function, x, hObj.ConfLevel);
    catch ignore %#ok<NASGU>
        % If errors occur in computing predictions and intervals, we
        % will still try to compute the predictions themselves below.
    end
end

% Compute predictions if the previous code didn't yield any
if isempty(y) && ~isempty(x)
    y = feval(hObj.Function, x);
end

% Replace any complex elements in the vectors to be plotted with NaNs.
y = i_NansFromComplexElements( y );
ci = i_NansFromComplexElements( ci );
end

function x = i_NansFromComplexElements( x )
% Replace each complex element, i.e., element with non-zero imaginary
% part, and set it to NaN.
%
% Determine which elements of the array have non-zero imaginary part.
tf = imag( x ) ~= 0;
% Set those elements to NaN
x(tf) = NaN;
end

function xlim = i_GetXLimFromAxesDataSpace( hObj )

ax = ancestor(hObj,'axes');
if isempty(ax)
    xlim = [];
    return;
end

hDataSpace = ax.DataSpace;
if isempty(hDataSpace)
    xlim = [];
    return;
end

xlim = hDataSpace.XLim;
end