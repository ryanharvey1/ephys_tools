function fitFigureDeleteFcn( this, ~, ~ )
%fitFigureDeleteFcn FitFigure DeleteFcn
%
%   fitFigureDeleteFcn(this, SOURCE, EVENT) is the FitFigure's DeleteFcn


%   Copyright 2008-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $    $Date: 2011/05/09 00:40:13 $

if isvalid(this.HSFTool)
    config = this.Configuration;
    config.Visible = 'off';
    this.HSFTool.HFitsManager.closeFit(this.FitUUID, config);
end
% delete the FittingPanel, so listeners do not try to respond
iDeleteListeners(this);
deletePanel(this.HFittingPanel);
delete( this.Handle );
delete( this );
% do the superclass close function
figureDeleteFcn(this);
end

function iDeleteListeners(this)
n = length(this.FitdevListeners);
for i = 1:n
    delete(this.FitdevListeners{i});
end
this.FitdevListeners = {};
end