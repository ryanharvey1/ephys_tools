function h = boundedline(fun, varargin)
% BOUNDEDLINE Create a functionline with confidence bounds.
%    H = BOUNDEDLINE(FUN) creates a line based on the function FUN. FUN
%    accepts vector input X and returns a vector output that is FUN evaluated
%    at X.
%
%    H = BOUNDEDLINE(FUN,SHOWBOUNDS,CONFLEVEL,DFE) controls confidence
%    bounds.  SHOWBOUNDS can be 'on' (the default) to show bounds or
%    'off' to omit them.  CONFLEVEL is the confidence level and is 0.95
%    by default.  DFE is the degrees of freedom for error and is Inf
%    by default.
%
%    H = BOUNDEDLINE(FUN,...,'-userargs',{P1,P2,...}) will pass the user 
%    arguments P1, P2,... as a comma separated list to FUN, e.g.
%    feval(FUN,X,P1,P2,...).
%
%    H = BOUNDEDLINE(FUN,...,'String','STR') specifies 'STR' as the
%    initial value of the String property.
%
%    H = BOUNDEDLINE(FUN,...,PROP1,VALUE1,PROP2,VALUE2,...)
%    will create a line with properties set to the given values,
%    e.g. 'color', 'b'.

%   Copyright 2000-2008 The MathWorks, Inc.
%   $Revision: 1.20.2.7 $  $Date: 2010/10/08 16:36:41 $

userargs = {};
showbounds = 'on';
conflevel = 0.95;
hstr = '';
dfe = Inf;
if nargin==0 || isequal(fun, 'parent')
   % This condition should only happen when we are being
   % called by the FIG-file loader which always passes
   % our constructor ('parent', parent) arguments, or
   % if we are otherwise being loaded from a file.
   fun = '';
   if nargin>1
      curraxes = varargin{1};
   else
      curraxes = gca;
   end
   varargin = {};
   userargs = {};
   loadFromFig = 1;
   cfitobj = cfit(fittype);
else
   loadFromFig = 0;
   dfe = fun.dfe;
   hstr = fun.name;
   cfitobj = fun.fit;
   
   if ~isempty(varargin)
      if isequal(varargin{1},'on') || isequal(varargin{1},'off')
         showbounds = varargin{1};
         varargin(1) = [];
      end
      if ~isempty(varargin) && isnumeric(varargin{1})
         if length(varargin{1})==1 && (varargin{1}>0) && (varargin{1}<1)
            conflevel = varargin{1};
            varargin(1) = [];
         else
            error(message('curvefit:boundedline:BadConfidenceLevel'));
         end
      end
      if ~isempty(varargin) && strcmpi(varargin{1},'string')
         hstr = varargin{2};
         varargin(1:2) = [];
      end
      if ~isempty(varargin) && strcmpi(varargin{1},'-userargs')
         userargs = varargin{2};
         varargin(1:2) = [];
      else
         userargs = {};
      end
   end
   
   parent = find(strncmpi('parent',varargin,6));
   if ~isempty(parent)
      curraxes = varargin{parent+1};
   else
      curraxes = gca;
   end
end

% Construct the function line (curve) - calls the built-in constructor
if isempty(varargin) && ~isempty(curraxes)
   h = cftool.boundedline('xdata',[], 'ydata',[], 'parent',curraxes);
else
   h = cftool.boundedline('xdata',[], 'ydata',[], varargin{:});
end

% Set FactoryValues here, since R12fcs ignores them
h.Granularity = 300;
h.UserArgs  = {};
h.String = '';
h.ShowBounds = 'on';
h.ConfLevel = 0.95;
h.DFE = Inf;
h.SSE = 0;
h.R = [];

% Initialize
h.function = cfitobj;   % Function to be evaluated
h.userargs = userargs;  % Cell array of user args
h.ShowBounds = showbounds;
h.ConfLevel = conflevel;
h.DFE = dfe;
h.String = hstr;
h.XLim = [];
h.fit = fun;
hh = [line(1,1,'Visible','off','Parent',curraxes), ...
      line(1,1,'Visible','off','Parent',curraxes)];
set(hh,'XData',[],'YData',[],'Visible','on','Tag','cfconf',...
       'LineStyle',':','UserData',h);
set(hh,'Color',get(h,'Color'),'ButtonDownFcn',get(h,'ButtonDownFcn'));
h.BoundLines = hh;

% Listeners
c = get(h,'Parent');
list = {
    % Listeners on this boundedline object, h
    curvefit.proplistener( h, 'UserArgs',    curvefit.callbackFunction( @localUpdateLine, h ) )
    curvefit.proplistener( h, 'Granularity', curvefit.callbackFunction( @localUpdateLine, h ) )
    curvefit.proplistener( h, 'String',      curvefit.callbackFunction( @localUpdateLine, h ) )
    curvefit.proplistener( h, 'ShowBounds',  curvefit.callbackFunction( @localUpdateLine, h ) )
    curvefit.proplistener( h, 'ConfLevel',   curvefit.callbackFunction( @localUpdateLine, h ) )
    curvefit.proplistener( h, 'DFE',         curvefit.callbackFunction( @localUpdateLine, h ) )
    % Listeners on the parent axes object
    curvefit.proplistener( c, 'XLim', curvefit.callbackFunction( @localUpdateLine, h ) )
    };
h.listeners = list;
set(h, 'DeleteFcn', curvefit.callbackFunction( @localDeleteBl, h ) );

if loadFromFig
   % This means we are loading from a file.
   % We want to turn off all of our listeners until we have been
   % fully constructed.  Then our createfunction which we
   % install here will re-enable the listeners.
   set( h, 'CreateFcn', curvefit.callbackFunction( @localCreateBl, h, list ) );
   for j=1:length(list)
      curvefit.setListenerEnabled( list{j}, false );
   end
else
    % Don't listen for axis limit changes during creation
    curvefit.setListenerEnabled( h.listeners{end}, false )
end

% Draw the line
update(h);

if ~loadFromFig
   % Start to listen now
   curvefit.setListenerEnabled( h.listeners{end}, true );
end

% Give it a context menu
c = findall(ancestor(curraxes,'figure'),'Type','uicontextmenu','Tag','fitcontext');
if ~isempty(c)
   set(double(h),'uiContextMenu',c);
end

%-----------------------------------------------------------
function localUpdateLine(src, ignore, hBline)
% LOCALUPDATELINE Callback to update boundedline when changes are required

% Watch out for spurious X limit notification
ax = get(hBline,'Parent');
chng = src.name;
bxlims = hBline.XLim;
if isequal(chng,'XLim') && ~isempty(bxlims) && isequal(get(ax,'XLim'),bxlims)
   return;
end

% Most changes require updating the curve, but string changes do not
if ~isequal(chng,'String')
   update(hBline);
end

% Some changes require updating the legend
if isequal(chng,'String') || isequal(chng,'ShowBounds')
   cffig = ancestor(hBline,'figure');
   
   % Update residual line so its legend will also update
   if isa(hBline.fit,'cftool.fit')
      rline = hBline.fit.rline;
      if ~isempty(rline) && ishandle(rline)
         ud = get(rline,'UserData');
         if iscell(ud)
            ud{1} = hBline.fit.name;
            set(rline,'UserData',ud);
         end
      end
   end
   cfswitchyard('cfupdatelegend',cffig);
end

%-----------------------------------------------------------
function localCreateBl(ignore1, ignore2, hBline, list)
for j=1:length(list)
      curvefit.setListenerEnabled( list{j}, true );
end
set(hBline, 'createfcn','');

%-----------------------------------------------------------
function localDeleteBl(ignore1, ignore2, hBline)
hh = hBline.BoundLines;
if ~isempty(hh)
   for j=1:length(hh)
      if ishandle(hh(j))
         delete(hh(j));
      end
   end
end
