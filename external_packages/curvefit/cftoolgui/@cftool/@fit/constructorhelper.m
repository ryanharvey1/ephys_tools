function h = constructorhelper(h,name)

% $Revision: 1.14.2.7 $  $Date: 2010/04/24 18:06:07 $
% Copyright 2001-2010 The MathWorks, Inc.

% This is the meat of the constructor.  It is here, rather than in fit.m,
% because there is no way to have another method, specifically copyfit, call the
% constructor.  Since it is a method of the fit object, it only calls the
% builtin.  An example of this workaround in action is in copyfit: 
% 
%   new = constructorhelper(cftool.fit);

NONE = cfswitchyard( 'cfGetNoneString' );

% Set the FactoryValues here, since R12fcs ignores them.
h.fit=cfit;
h.isGood=false;
h.fitOptions=curvefit.basefitoptions;
h.outlier= NONE;
h.goodness=[];
h.sse=NaN;
h.rsquare=NaN;
h.adjsquare=NaN;
h.rmse=NaN;
h.dfe=NaN;
h.ncoeff=NaN;
h.output=[];
h.plot=0;
h.line=[];
h.rline=[];
h.ColorMarkerLine = [];

if nargin==1,
    taken = 1;
    count=cfgetset('fitcount');
    if isempty(count)
        count = 1;
    end
    while taken
        name=sprintf('fit %i', count);
        if isempty(find(getfitdb,'name',name))
            taken = 0;
        else
            count=count+1;
        end
    end
    cfgetset('fitcount',count+1);
end

h.name=name;

list(1) = handle.listener(h,findprop(h,'fit'),'PropertyPostSet', {@fitplot,h});
list(2) = handle.listener(h,findprop(h,'plot'),'PropertyPostSet', {@fitplot,h});
list(3) = handle.listener(h,'ObjectBeingDestroyed', {@cleanup,h});
list(4) = handle.listener(h,findprop(h,'dshandle'),'PropertyPostSet',...
                          {@updatelim,h});
list(5) = handle.listener(h,findprop(h,'name'),'PropertyPostSet', {@newname,h});

h.listeners=list;

cfgetset('dirty',true);   % session has changed since last save

%=============================================================================
function updatelim(~,~,fit)
% updatelim( src, evt, fit )
fit.xlim = xlim(fit);

%=============================================================================
function newname(~, ~,fit)
% newname(src, evt ,fit)
h = fit.line;
if isa(h,'cftool.boundedline') || isa(h,'cftool.BoundedFitLine') 
   h.String = fit.name;
end
cfgetset('dirty',true);   % session has changed since last save

%=============================================================================
function fitplot(~, ~,fit)
% fitplot( src, evt, fit )
updateplot(fit);
cfgetset('dirty',true);   % session has changed since last save

%=============================================================================
function cleanup(~, ~,fit)
% cleanup(src, evt, ,fit)
if ~isempty(fit.line) && ishandle(fit.line)
   fit.plot = 0;
end
if ~isempty(fit.rline) && ishandle(fit.rline)
   delete(fit.rline);
end
cfgetset('dirty',true);   % session has changed since last save
