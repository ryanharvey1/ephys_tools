function setLimits( this, ax, limitPVPairs )
% setLimits sets the AxesViewModel Limit properties
%
%   setLimits(this, ax, limitPVPairs) sets the AxesViewModel Limit
%   properties based on the limits of AX. limitPVPairs is a cell array of
%   axes Limit properties and their values.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.4.1 $    $Date: 2011/07/18 00:31:16 $

[xlim, ylim, zlim] = iGetLimits(limitPVPairs{:});

if iIsResidualsAxes(ax)
    if isCurveDataSpecified( this.HFitdev.FittingData )
        setLimits(this.AxesViewModel, {xlim}, [], ylim);
    else
        setLimits(this.AxesViewModel, {xlim, ylim}, [], zlim);
    end
else % main or contour axes
    if isCurveDataSpecified( this.HFitdev.FittingData )
        setLimits(this.AxesViewModel, {xlim}, ylim, []);
    else
        setLimits(this.AxesViewModel, {xlim, ylim}, zlim, []);
    end
end
end

function [xlim, ylim, zlim] = iGetLimits(varargin)
% iGetLimits extracts limits from the name/value pairs of varargin
xlim = [];
ylim = [];
zlim = [];
for i = 1:2:length( varargin )
    switch lower(varargin{i})
        case 'xlim'
            xlim = varargin{i+1};
        case 'ylim'
            ylim = varargin{i+1};
        case 'zlim'
            zlim = varargin{i+1};
        otherwise
            warning(message('curvefit:sftoolgui:FitFigure:setLimits:UnexpectedLimit', varargin{i}));
    end
end
end

function tf = iIsResidualsAxes(ax)
% iIsResidualsAxes returns true if the axes is from the residuals plot and
% false otherwise.
tf = strcmpi(get(ax, 'Tag'), 'sftool residuals axes');
end

