function [emsg,Xstr,Lstr,Fstr,Ustr,Dstr,D2str,Istr,wmsg] = cfanalysis(expr,analysis)
% CFANALYSIS Helper function for the Curve Fitting toolbox Analysis panel
%
%    [XSTR,LSTR,FSTR,USTR,DSTR,D2STR,ISTR] = CFANALYSIS(EXPRESSION,
%    ANALYSIS)
%    evaluates EXPRESSION to produce values of the independent
%    variable X, queries the ANALYSIS object to determine what
%    calculations are required, and returns the required results.
%    The outputs are string arrays for X, lower bound on fit, fit,
%    upper bound on fit, 1st derivative, 2nd derivative, integral.

%   Copyright 2001-2009 The MathWorks, Inc.
%   $Revision: 1.24.2.9 $  $Date: 2009/05/07 18:17:42 $

Xstr = '';
Lstr = '';
Fstr = '';
Ustr = '';
Dstr = '';
D2str = '';
Istr = '';
wmsg = '';

[lw,lwid] = lastwarn;
ws = warning('off', 'all');
lastwarn('');

try
    [emsg,Xstr,Lstr,Fstr,Ustr,Dstr,D2str,Istr,wmsg] = cfsub(expr,analysis);
catch e
    warning(ws);
    emsg = sprintf('Error:\n%s',e.message);
end

lastwarn(lw,lwid);
warning(ws);

%--------------------------------------------------------------------------
function [emsg,Xstr,Lstr,Fstr,Ustr,Dstr,D2str,Istr,wmsg] = cfsub(expr,analysis)
emsg = '';
Xstr = '';
Lstr = '';
Fstr = '';
Ustr = '';
Dstr = '';
D2str = '';
Istr = '';
wmsg = '';

% Clear plot if requested
plotfig = cfgetset('analysisfigure');
if isequal(analysis,'clear')
    if ~isempty(plotfig) && ishandle(plotfig)
        h = findobj(allchild(plotfig),'flat','serializable','on');
        delete(h);
    end
    h = axes('Visible','off','Parent',plotfig);
    text(.5,.5,xlate('Press "Apply" to create a new plot'),...
        'Parent',h,'HorizontalAlignment','center');
    return
end

% Find the fit used for the evaluation, get its cfit
fitobj = handle(analysis.getFit);
if isempty(fitobj)
    emsg = sprintf('Cannot find fit.');
    emsg = combinewe(wmsg,emsg);
    return;
end
ds = fitobj.dshandle;

% Get a set of x values at which to evaluate stuff
inclvec = ~getexcluded(ds,fitobj.outlier);
xdata = ds.x(inclvec);
ydata = ds.y(inclvec);
% waitbar is shown only if array lengths are greater than 15000
hWaitbar = [];
try
    x = cfeval(xdata, ydata, [], expr);
    wmsg = combineww(wmsg);
    if ~isempty(x) && length(x) > 15000  %put up a waitbar
        hWaitbar = waitbar(0,'Please wait...', 'Name', 'Data Loading');
        waitbar(.05, hWaitbar);
    end
catch e
    emsg = sprintf('Attempt to evaluate Xi expression failed with the following error:\n\n%s',...
        e.message);
    emsg = combinewe(wmsg,emsg);
    return;
end
if isempty(x)
    emsg = sprintf('Xi expression yielded no data.');
    emsg = combinewe(wmsg,emsg);
    return;
end
x = sort(x(:));
if ~isempty(hWaitbar) && ishandle(hWaitbar)
    waitbar(.10, hWaitbar);
end

tmp = sprintf('%10g\n',x); 
Xstr = strread(tmp,'%s');
if ~isempty(hWaitbar) && ishandle(hWaitbar)
    waitbar(.15, hWaitbar);
end
results.xi = x;
if ~isempty(hWaitbar) && ishandle(hWaitbar)
    waitbar(.20, hWaitbar);
end
% Find the fit used for the evaluation, get its cfit
cfitobj = get(fitobj,'fit');

% If requested, evaluate the fit and compute bounds for it
yfit = [];
ybnds = [];

% Default title if no curve plotted
t1 = 'Data';
if ~isempty(hWaitbar) && ishandle(hWaitbar)
    waitbar(.25, hWaitbar);
end
level = 95;
if analysis.isInterpolateSelected && ~isempty(x)
    % First get confidence level if needed
    bounds = analysis.getInterpolateBounds;
    if bounds~=0
        level = char(analysis.getInterpolateLevel);
        try
            level = cfeval(xdata, ydata, x, level);
            wmsg = combineww(wmsg);
        catch e
            emsg = sprintf('Attempt to evaluate Level expression failed with the following error:\n\n%s',...
                e.message);
            emsg = combinewe(wmsg,emsg);
            if ~isempty(hWaitbar) && ishandle(hWaitbar)
                close(hWaitbar)
            end
            return;
        end
        if ~isnumeric(level) || length(level)~=1 || level<=0 || level>=100
            emsg = 'Confidence level must be positive and less than 100.';
            emsg = combinewe(wmsg,emsg);
            if ~isempty(hWaitbar) && ishandle(hWaitbar)
                close(hWaitbar)
            end
            return;
        end
    end
    
    if ~isempty(hWaitbar) && ishandle(hWaitbar)
        waitbar(.30, hWaitbar);
    end

    % Get predictions and bounds
    try
        switch bounds
            case 0
                yfit = feval(cfitobj,x);
                t1 = 'Fit';
            case 1
                [ybnds,yfit] = predint(cfitobj,x,level/100,'functional');
                t1 = sprintf('Fit with %g%% pred bounds',level);
                results.boundtype = 'functional';
            case 2
                [ybnds,yfit] = predint(cfitobj,x,level/100,'observation');
                t1 = sprintf('Fit with %g%% pred bounds',level);
                results.boundtype = 'observation';
        end
        wmsg = combineww(wmsg);
    catch e
        emsg = sprintf(...
            'Error evaluating fit.  Error message is:\n%s',...
            e.message);
        emsg = combinewe(wmsg,emsg);
        if ~isempty(hWaitbar) && ishandle(hWaitbar)
            close(hWaitbar)
        end
        return;
    end
    if ~isempty(ybnds)
        if ~isempty(hWaitbar) && ishandle(hWaitbar)
            waitbar(.35, hWaitbar);
        end
        tmp = sprintf('%6g\n',ybnds(:,1)); 
        Lstr = strread(tmp,'%s');
        tmp = sprintf('%6g\n',ybnds(:,2)); 
        Ustr = strread(tmp,'%s');
        if ~isempty(hWaitbar) && ishandle(hWaitbar)
            waitbar(.40, hWaitbar);
        end
        results.lower = ybnds(:,1);
        results.upper = ybnds(:,2);
        results.conflevel = level/100;
    end
    if ~isempty(yfit)
        if ~isempty(hWaitbar) && ishandle(hWaitbar)
            waitbar(.45, hWaitbar);
        end
        tmp = sprintf('%6g\n',yfit); 
        Fstr = strread(tmp,'%s');
        if ~isempty(hWaitbar) && ishandle(hWaitbar)
            waitbar(.50, hWaitbar);
        end
        results.yfit = yfit;
    end
end
if ~isempty(hWaitbar) && ishandle(hWaitbar)
     waitbar(.55, hWaitbar);
end
% If requested, evaluate the derivatives
yderiv = [];
yderiv2 = [];
try
    if analysis.isSecondDerivativeSelected
        [yderiv,yderiv2] = differentiate(cfitobj,x);
        results.d2ydx2 = yderiv2;
    end
    if analysis.isFirstDerivativeSelected
        if isempty(yderiv)
            yderiv = differentiate(cfitobj,x);
        end
        results.dydx = yderiv;
    else
        yderiv = [];
    end
    wmsg = combineww(wmsg);
catch e
    emsg = sprintf(...
        'Error computing derivative.  Error message is:\n%s',...
        e.message);
    emsg = combinewe(wmsg,emsg);
    if ~isempty(hWaitbar) && ishandle(hWaitbar)
        close(hWaitbar)
    end
    return;
end
if ~isempty(yderiv)
    if ~isempty(hWaitbar) && ishandle(hWaitbar)
        waitbar(.60, hWaitbar);
    end
    tmp = sprintf('%6g\n',yderiv); 
    Dstr = strread(tmp,'%s');
end
if ~isempty(yderiv2)
    if ~isempty(hWaitbar) && ishandle(hWaitbar)
         waitbar(.65, hWaitbar);
    end
    tmp = sprintf('%6g\n',yderiv2); 
    D2str = strread(tmp,'%s');
end
if ~isempty(hWaitbar) && ishandle(hWaitbar)
    waitbar(.70, hWaitbar);
end

% If requested, evaluate the integral
yint = [];
if analysis.isIntegrateSelected
    if analysis.isStartMinXSelected
        x0 = min(x);
    else
        start = char(analysis.getIntegrateStartValue);
        try
            x0 = cfeval(xdata, ydata, x, start);
            wmsg = combineww(wmsg);
        catch e
            emsg = sprintf('Attempt to evaluate "Start from" expression failed with the following error:\n\n%s',...
                e.message);
            emsg = combinewe(wmsg,emsg);
            if ~isempty(hWaitbar) && ishandle(hWaitbar)
                close(hWaitbar)
            end
            return;
        end
        if ~isnumeric(x0) || length(x0)~=1 || isnan(x0) || ~isreal(x0)
            emsg = 'Starting point must be a real number.';
            emsg = combinewe(wmsg,emsg);
            if ~isempty(hWaitbar) && ishandle(hWaitbar)
                close(hWaitbar)
            end
            return;
        end
    end
    if ~isempty(hWaitbar) && ishandle(hWaitbar)
        waitbar(.75, hWaitbar);
    end

    try
        yint = integrate(cfitobj,x,x0);
        t4 = sprintf('Integral from %g',x0);
        results.integral = yint;
        results.integralstart = x0;
        wmsg = combineww(wmsg);
    catch e
        emsg = sprintf(...
            'Error computing integral.  Error message is:\n%s',...
            e.message);
        emsg = combinewe(wmsg,emsg);
        if ~isempty(hWaitbar) && ishandle(hWaitbar)
            close(hWaitbar) 
        end
        return;
    end
    if ~isempty(yint)
        if ~isempty(hWaitbar) && ishandle(hWaitbar)
             waitbar(.80, hWaitbar);
        end
        tmp = sprintf('%6g\n',yint); 
        Istr = strread(tmp,'%s');
    end
end
if ~isempty(hWaitbar) && ishandle(hWaitbar)
    waitbar(.85, hWaitbar);
end
% Now deal with plotting
plotfig = cfgetset('analysisfigure');

plotpoints = analysis.isPlotDataSetSelected;
plotresults = analysis.isPlotResultsSelected;
nplots = (plotpoints | (plotresults & ~isempty(yfit))) + ...
    plotresults * (~isempty(yderiv) + ~isempty(yderiv2) + ~isempty(yint));
if nplots > 0
    % First see if there is a desired color
    fcolor = [0 0 1];   % blue
    dcolor = [1 0 0];   % red
    h = fitobj.line;
    if ~isempty(h) && ishandle(h)
        fcolor = get(h,'Color');
    end
    ds = fitobj.dshandle;
    if ~isempty(ds) && ishandle(ds)
        h = ds.line;
        if ~isempty(h) && ishandle(h)
            dcolor = get(h,'Color');
        end
    end

    % Create plotting figure if it does not yet exist
    if isempty(plotfig) || ~ishandle(plotfig)
        plotfig=figure(...
            'IntegerHandle','off',...
            'HandleVisibility','callback',...
            'name','Curve Fitting Analysis',...
            'numbertitle','off',...
            'PaperPositionMode','auto');
        cfgetset('analysisfigure',plotfig);
    end

    % New or old, prepare figure by removing old contents
    h = findobj(allchild(plotfig),'flat','serializable','on');
    delete(h);
    if ~isempty(hWaitbar) && ishandle(hWaitbar)
        waitbar(.95, hWaitbar);
    end
    % Do each requested plot
    k = 1;
    allax = zeros(nplots,1);
    if plotpoints || ~isempty(yfit)
        ax = localsubplot(nplots,k,plotfig);
        allax(k) = ax;
        h = zeros(3,1);
        if plotresults && ~isempty(yfit)
            h(1) = plot(x,yfit,'.-','Parent',ax,'Color',fcolor);
        end
        if plotresults && ~isempty(ybnds)
            h(2)=line(x,ybnds(:,1),'Parent',ax,'Color',fcolor,'LineStyle',':');
            line(x,ybnds(:,2),'Parent',ax,'Color',fcolor,'LineStyle',':');
        end
        if plotpoints
            h(3) = line(xdata,ydata,'Parent',ax,...
                'Color',dcolor,'LineStyle','none','Marker','x');
        end
        setylabel(ax,t1);
        k = k+1;

        % Add a legend to this plot
        if any(h~=0)
            c = cell(3,1);
            c{1} = fitobj.name;
            c{2} = sprintf('%g%% prediction bounds',level);
            c{3} = ds.name;
            c(h==0) = [];
            h(h==0) = [];
            legh = legend(h,c);
            set(legh,'Interpreter','none');
        end
    end
    if plotresults && ~isempty(yderiv)
        ax = localsubplot(nplots,k,plotfig);
        allax(k) = ax;
        plot(x,yderiv,'.-','Parent',ax,'Color',fcolor);
        setylabel(ax,'1st deriv');
        k = k+1;
    end
    if plotresults && ~isempty(yderiv2)
        ax = localsubplot(nplots,k,plotfig);
        allax(k) = ax;
        plot(x,yderiv2,'.-','Parent',ax,'Color',fcolor);
        setylabel(ax,'2nd deriv');
        k = k+1;
    end
    if plotresults && ~isempty(yint)
        ax = localsubplot(nplots,k,plotfig);
        allax(k) = ax;
        plot(x,yint,'.-','Parent',ax,'Color',fcolor);
        setylabel(ax,t4);
    end

    % Tick labels on upper axes just get in the way
    if nplots>1
        set(allax(1:end-1),'XTickLabel','');
        xlim = get(allax,'XLim');
        xlim = cat(2,xlim{:});
        xlim = [min(xlim) max(xlim)];
        set(allax,'XLim',xlim);
    end

    % Add title to top graph
    txt = sprintf('Analysis of fit "%s" for dataset "%s"',...
        fitobj.name, ds.name);
    h = get(allax(1),'Title');
    set(h,'String',txt,'Interpreter','none');

    figure(plotfig);

    % If no plotting selected
elseif ~isempty(plotfig) && ishandle(plotfig)
    delete(plotfig);
end
if ~isempty(hWaitbar) && ishandle(hWaitbar)
    waitbar(1, hWaitbar);
end
% Save results so they can be saved to workspace if requested later
cfgetset('analysisresults',results);

% Clean up
wmsg = combineww(wmsg);

if ~isempty(hWaitbar) && ishandle(hWaitbar)
    close(hWaitbar)
end


% -------------  Set up axes to plot kth out of n plots
function ax = localsubplot(n,k,f)

pos = [0.13 0.11 0.775 0.815];
h = pos(4) / (n + (n-1)/10);     % height of each plot
g = h/10;                        % gutter between plots
pos(2) = pos(2) + pos(4) - k*h - (k-1)*g;
pos(4) = h;
ax = axes('Position',pos,'Parent',f);


% -------------- Set y axis label
function setylabel(ax,txt)

h = get(ax,'YLabel');
set(h,'String',txt);

% -------------- Evaluate text field in defined environment
function varargout = cfeval(x,y,xi,varargin) %#ok use x, y, and xi if they
                                             % are part of the expression
                                             % in varargin
varargout{1} = eval(['[' varargin{1} ']']);


% -------------- Combine warning and error messages
function e = combinewe(w,e)

if ~isempty(w)
    e = sprintf('%s\n\nThe following warning also occurred:\n%s',xlate(e),xlate(w));
end

% -------------- Combine existing warning msg with current warning, clear current
function w = combineww(w)

if ~isempty(lastwarn)
    if isempty(w)
        w = sprintf('Warning:\n%s',lastwarn);
    elseif isempty(strfind(w,lastwarn))
        w = sprintf('%s\n\n%s',w,lastwarn);
    end
    lastwarn('');
end
