function cfupdatelegend(cffig,reset,leginfo,rleginfo)
%CFUPLDATELEGEND Update legend in curve fitting plot

%   $Revision: 1.20.2.11 $  $Date: 2010/05/10 16:59:42 $ 
%   Copyright 2000-2010 The MathWorks, Inc.

if nargin<2, reset=false; end
if nargin<3, leginfo={}; end
if nargin<4, rleginfo={}; end

% If figure not passed in, find figure that contains this thing
cffig = ancestor(cffig,'figure');

% Remember info about old legend, if any
ax = findall(cffig,'Type','axes','Tag','main');
[leginfo,havelegendloc,relpos] = localgetlegendinfo(cffig,ax,reset,leginfo);
legend(ax, 'off');

% Remember info about old residuals legend, if any
ax2 = findall(cffig,'Type','axes','Tag','resid');
if ~isempty(ax2)
   [rleginfo,haverlegendloc,rrelpos] = localgetlegendinfo(cffig,ax2,reset,rleginfo);
   legend(ax2, 'off');
end

% Maybe no legend has been requested
if isequal(cfgetset('showlegend'),'off')
   return
end

% Get data line handles and labels
hh = flipud(findobj(ax,'Type','line'));
h1 = findobj(hh,'flat','Tag','cfdata');
n = length(h1);
c1 = cell(n,1);
toKeep = true( size( h1 ) );
for j=1:length(h1)
   nm = '';
   ud = get(h1(j),'UserData');
   if ~isempty(ud) && ishandle(ud) && ~isempty(findprop(ud,'name'))
      nm = ud.name;
   end
   if isempty(nm)
       toKeep(j) = false;
   else
      c1{j} = nm;
   end
end
c1 = c1(toKeep);
h1 = h1(toKeep);
s1 = 1000*(1:length(h1));
if isempty(s1)
   maxnum = 0;
else
   maxnum = max(s1) + 1000;
end

% Indent fits if there are two or more data lines
if (length(h1)>1)
   pre = '  ';
else
   pre = '';
end

% Get fit line handles and labels
h2 = findobj(hh,'flat','Tag','curvefit');
s2 = zeros( size( h2 ) ); % index of position in legend
n = length(h2);
c2 = cell(n,1);
nms = cell(n,1);
toKeep = true( size( h2 ) );
for j=1:length(h2)
   nm = get(handle(h2(j)),'String');
   if isempty(nm)
       toKeep(j) = false;
   else
      nms{j} = nm;
      c2{j} = [pre nm];
      
      % Find the dataset for this fit
      ua = get(handle(h2(j)),'UserArgs');
      if iscell(ua) && ~isempty( ua )
         ds = ua{1};
         s2j = maxnum + j;
         for k=1:length(h1)
            if isequal(ds,c1{k})
               s2j = s1(k) + j;
               break;
            end
         end
         s2(j) = s2j;
      end
   end
end
nms = nms(toKeep);
c2 = c2(toKeep);
h2 = h2(toKeep);
s2 = s2(toKeep);

% Remember just one confidence bound for each fit
hconf = zeros(size(h2));
for j=1:length(hconf)
   b = get(handle(h2(j)),'BoundLines');
   if ~isempty(b)
      hconf(j) = b(1);
   end
end

% Indent bounds if there are two or more fits
if (length(h2)>1)
   pre = [pre '  '];
end

% Get confidence bound line handles and labels
h3 = findobj(hh,'flat','Tag','cfconf');
n = length(h3);
c3 = cell(n,1);
s3 = zeros(size(h3));
toKeep = true( size( h3 ) );
for j=1:length(h3)
   if isempty(hconf)
      k = [];
   else
      k = find(h3(j)==hconf);
   end
   if ~isempty(k) && ~isempty(get(h3(j),'XData'))
      c3{j} = sprintf('%sPred bnds (%s)',pre,nms{k});
      s3(j) = s2(k) + 0.5;
   else
       toKeep(j) = false;
   end
end
c3 = c3(toKeep);
h3 = h3(toKeep);
s3 = s3(toKeep);

% Combine everything together for the legend
h = [h1(:); h2(:); h3(:)];
c = [c1; c2; c3];
s = [s1(:); s2(:); s3(:)];

% Sort so related things are together
[~,j] = sort(s);
c = c(j);
h = h(j);

% Create the legend
if ~isempty( h )
   ws = warning;
   lw = lastwarn;
   warning('off', 'all');
   legh = [];
   try
      legh = legend(ax,h,c,leginfo{:});
      if ~havelegendloc && ~isempty(relpos)
          setrelativelegendposition(relpos,cffig,ax,legh);
      end
      localFixContextMenu( legh );
   catch ignore %#ok<NASGU>
   end
   warning(ws);
   lastwarn(lw);
   
   % Avoid treating ds/fit names as TeX strings
   set(legh,'Interpreter','none');
else
   legend(ax,'off');
end

% Set a resize function that will handle legend and layout
set(cffig,'resizefcn','cftool(''adjustlayout'');');

% Fix residual legend; this is a lot simpler
if ~isempty(ax2)
   h = flipud(findobj(ax2,'Type','line'));
   c = cell(length(h),1);
   for j=length(h):-1:1
      t = get(h(j),'UserData');
      if iscell(t) && ~isempty( t ) && ischar(t{1}) && ~isempty(t{1})
         c{j} = t{1};
      else
         c(j) = [];
         h(j) = [];
      end
   end
   
   if ~isempty( h )
      ws = warning;
      lw = lastwarn;
      warning('off', 'all');
      legh = [];
      try
         legh = legend(ax2,h,c,rleginfo{:});
         if ~haverlegendloc && ~isempty(rrelpos)
            setrelativelegendposition(rrelpos,cffig,ax2,legh);
         end
         localFixContextMenu( legh );
      catch ignore %#ok<NASGU>
      end
      warning(ws);
      lastwarn(lw);

      % Avoid treating ds/fit names as TeX strings
      if ~isempty(legh)
         set(legh,'Interpreter','none');
      end
   else
      legend(ax2,'off');
   end
end

% ---------------------------------------------------------
function [leginfo,havelegendloc,relpos] = localgetlegendinfo(cffig,ax,reset,leginfo)
% Get legend information

legh = legend(ax);
relpos = [];
if isempty(leginfo)
    if ~isempty(legh) && ishandle(legh) && ~reset
        leginfo = cfgetlegendinfo(legh);
    else
        leginfo = {};
    end
end

% Loop to find 'Location', as some non-text entries make ismember fail
havelegendloc = false;
for j=1:length(leginfo)
    if isequal('Location',leginfo{j})
        havelegendloc = true;
        break
    end
end
if ~havelegendloc && ~isempty(legh) && ishandle(legh) && ~reset
   relpos = getrelativelegendposition(cffig,ax,legh);
end

% ---------------------------------------------------------
function localFixContextMenu( hLegend )
% The legend gets created with a context menu. However this context menu
% has some features that have a destructive affect on CFTOOL. In this
% little function, we remove those features....
cmh = get( hLegend, 'UIContextMenu' );
% The children (menu entries) of the context menu are hidden so we need
% to get around that
h = allchild( cmh );
% Our actions are based on tags of items that appear in the context menu so
% we need to get all of those tags.
tags = get( h, 'Tag' );

% Delete the entries that cause bad things to happen
% Find tags of entries to delete.
TAGS_TO_DELETE = {'scribe:legend:mcode', 'scribe:legend:propedit', 'scribe:legend:interpreter'};
tf = ismember( tags, TAGS_TO_DELETE );
delete( h(tf) );

% For the 'Delete' item, we want to redirect the call to the CFTOOL
% legend toggle function
tf = ismember( tags, 'scribe:legend:delete' );
set( h(tf), 'Callback', @(s, e) cftool( 'togglelegend', 'off' ) );

% For the 'Refresh' item, we want to redirect the callback to reset the
% properties of the legend, e.g., the colour and font.
tf = ismember( tags, 'scribe:legend:refresh' );
set( h(tf), 'Callback', @(s, e) cfupdatelegend( cfgetset( 'cffig' ), true ) );

