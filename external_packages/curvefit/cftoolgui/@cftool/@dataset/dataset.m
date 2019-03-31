function h = dataset(x,y,weight,dsname,xvals,yvals,wvals,src)
%DATASET Constructor for CFTOOL dataset
%
%   H = DATASET( 'disconnected' ) creates dataset that is not connected to the
%   dataset database.
%
%   H = DATASET(X, Y, WEIGHT, DSNAME, XVALS, YVALS, WVALS, SRC)

% $Revision 1.34.6.4 $  $Date: 2009/12/02 06:41:39 $
% Copyright 2000-2009 The MathWorks, Inc.

NONE= cfswitchyard( 'cfGetNoneString' );

h = cftool.dataset;

% Define FactoryValues here, since R12fcs ignores them
h.xname=NONE;
h.yname=NONE;
h.weightname=NONE;
h.plot=1;
h.ColorMarkerLine = [];


nargs = nargin;

% We may be asked to create an empty object not connected to the database
toconnect = 1;
if nargs==1 && isequal(x,'disconnected')
   toconnect = 0;
   h.plot=0;
   nargs = 0;
end

% For constructing empty object
if nargs==0
   x = 'x';
   y = 'y';
   weight = NONE;
   dsname = '';
   xvals = [];
   yvals = [];
   wvals = [];
   src = [];
   nargs = 8;
end

% if only x is given, swap x and y
if isequal(y,NONE)
   y=x;
   if nargin>=5
      yvals = xvals;
      nargs = max(6,nargs);
   end
   x=NONE;
end

if nargs>=6
   h.y = yvals;
else
   h.y=evalin('base',y);
end
if isempty(y)
    y = 'y';
end

if isequal(x,NONE)
   h.x=1:length(h.y);
elseif nargs>=5
   h.x = xvals;
else
   h.x=evalin('base',x);
end
if isempty(x)
    x = 'x';
end

if isequal(weight, NONE)
	h.weight=[];
elseif nargs>=7
   h.weight = wvals;
else
	h.weight=evalin('base',weight);
end
if isempty(weight)
    weight = 'w';
end


if nargs<4 || isempty(dsname)
   if isequal(x,NONE)
      dsname=y;
   else
      dsname=sprintf('%s vs. %s',y,x);
   end
   if ~isequal(weight,NONE)
      dsname = sprintf('%s with %s',dsname,weight);
   end
end

%make sure default name is unique 
if ~isempty(find(getdsdb,'name',dsname))
   taken = 1; 
   i = 2; 
   % search for first unique name 
   while taken 
      tryname = sprintf('%s (%i)', dsname, i); 
      if isempty(find(getdsdb,'name',tryname))
         dsname = tryname;
         taken = 0; 
      else 
         i=i+1; 
      end 
   end 
end
h.name = dsname;

if nargs>=8
   h.source = src;
else
   h.source=[];
end


% Make sure we store data as column vectors
h.x = h.x(:);
h.y = h.y(:);
h.weight = h.weight(:);
h.xlength = length(h.x);

% Now check for NaN, Inf or complex values in h.x, h.y and h.weight
xx=h.x;
yy=h.y;
ww=h.weight;

nans = false;
if any(isnan(xx))
    nans = true;
end
if any(isnan(yy))
    nans = true;
end
if any(isnan(ww))
    nans = true;
end

cplx = false;
if ~isreal(xx) 
    xx = real(xx);
    cplx = true;
end
if ~isreal(yy)
    yy = real(yy);
    cplx = true;
end
if ~isreal(ww)
    ww = real(ww);
    cplx = true;
end

infs = false;
if any(isinf(xx))
    xx(isinf(xx)) = NaN;
    infs = true;
end
if any(isinf(yy))
    yy(isinf(yy)) = NaN;
    infs = true;
end
if any(isinf(ww))
    ww(isinf(ww)) = NaN;
    infs = true;
end

h.x = xx;
h.y = yy;
h.weight=ww;

import com.mathworks.toolbox.curvefit.Analysis;

% using Analysis as parent to get the MATLAB icon. It is always initialized.
% This is required when starting cftool with data. 

errstring = '';
if nans
    errstring = xlate( 'Ignoring NaNs in data.' );
end
if infs
    inferr = xlate( 'Ignoring Infs in data.' );
    if isempty(errstring)
        errstring = inferr;
    else
        errstring = sprintf('%s\n%s', errstring, inferr);
    end
end
if cplx
    cplxerr = xlate( 'Using only the real component of complex data.' );
    if isempty(errstring)
        errstring = cplxerr;
    else
        errstring = sprintf('%s\n%s', errstring, cplxerr);
    end
end

if ~isempty(errstring)
     uiwait(msgbox(errstring, 'Import Data', 'warn', 'modal'));
end

% Remember the names
h.xname = x;
h.yname = y;
h.weightname = weight;

updatelim(h);
% add listeners
list(1) = handle.listener(h,findprop(h,'plot'),'PropertyPostSet', {@update,h});
list(2) = handle.listener(h,'ObjectBeingDestroyed', {@cleanup,h});
list(3) = handle.listener(h,findprop(h,'name'),'PropertyPostSet', {@newname});
h.listeners=list;

% add it to the list of datasets
if toconnect
   connect(h,getdsdb,'up');
   cfgetset('dirty',true);   % session has changed since last save
end

if h.plot==1 && nargin>0
    updatelim(h);
    cfswitchyard('cfmplot',h);
end

% ---- listener updates limits and calls plotting function
function update(~,~,ds)
% UPDATE(SRC, EVT, DS)
updatelim(ds);
dsmgr = com.mathworks.toolbox.curvefit.DataSetsManager.getDataSetsManager;
dsmgr.dataSetListenerTrigger(java(ds), dsmgr.DATA_SET_CHANGED, '', '');
cfswitchyard('cfmplot',ds);
cfgetset('dirty',true);   % session has changed since last save


% ---- listener updates name in legend
function newname(~,~)
% NEWNAME(SRC, EVT)
cfswitchyard('cfupdatelegend',cfgetset('cffig'));
cfgetset('dirty',true);   % session has changed since last save

% ---- Remember x and y limits, useful for selecting plot limits
function updatelim(h)

if isempty(h.x)
   h.xlim = [];
else
   h.xlim = [min(h.x) max(h.x)];
end
if isempty(h.y)
   h.ylim = [];
else
   h.ylim = [min(h.y) max(h.y)];
end

% ---- listener to unplot dataset when the object is destroyed
function cleanup(~,~,dataset)
% CLEANUP(SRC,EVT,DATASET)
changed = 0;
if ~isempty(dataset.line)
   if ishandle(dataset.line)
      list = dataset.listeners;
      list(1).enable = 'off'; %#ok<NASGU> handle class so "list" is used.
      dataset.plot = 0;
      changed = 1;
      cfswitchyard('cfmplot',dataset);
   end
end
if changed
   cfswitchyard('cfupdatelegend',cfgetset('cffig'));
   cfgetset('dirty',true);   % session has changed since last save
end
