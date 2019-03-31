function cfmplot(ds)
%CFMPLOT Create plot for curve fitting gui

%   $Revision: 1.14.2.6 $  $Date: 2007/06/14 04:55:04 $
%   Copyright 2000-2007 The MathWorks, Inc.


% Re-create figure if it does not exist
cffig = cftool('makefigure');

ax = findall(cffig,'Type','axes','Tag','main');

if ds.plot
    if isempty(ds.line) || ~ishandle(ds.line)
        if isequal(size(ds.x(:)), size(ds.y(:)))
           x = ds.x;
           y = ds.y;
        else
           x = [];
           y = [];
        end
        [c,m,l,w] = cfswitchyard('cfgetcolor',ax,'data',ds);
        ds.line=line('XData',x, 'YData',y, 'Parent',ax,'ButtonDownFcn',@cftips,...
                     'marker',m,'linestyle',l,'color',c,...
                     'LineWidth',w,'tag','cfdata','UserData',ds);
        
        if isequal(m,'.')
           set(ds.line, 'MarkerSize',12);
        end
        savelineproperties(ds);

        % Give it a context menu
        c = findall(ancestor(ax,'figure'),'Type','uicontextmenu',...
                    'Tag','datacontext');
        if ~isempty(c)
           set(ds.line,'uiContextMenu',c);
        end

        if ~isempty(ds.x)
           cfupdatexlim([min(ds.x) max(ds.x)]);
           cfupdateylim;
        end
        cfupdatelegend(cffig);
    end
else
    if (~isempty(ds.line)) && ishandle(ds.line)
        savelineproperties(ds);
        
        delete(ds.line);
        ds.line=[];
        zoom(cffig,'reset');
        
        % Update axis limits
        if ~isempty(ds.x)
           cfupdatexlim;
           cfupdateylim;
        end

        cfupdatelegend(cffig);
    end
end
