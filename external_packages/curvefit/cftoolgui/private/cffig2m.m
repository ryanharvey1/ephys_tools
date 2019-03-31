function cffig2m( varargin )
%CFFIG2M  Generate MATLAB code from the current CFTOOL session
%
%   CFFIG2M and CFFIG2M(CFFIG) open a new document in the MATLAB editor and
%   write code to that document.
%
%   CFFIG2M(CFFIG, FILENAME) creates a new file, FILENAME, and write code to
%   that file.

%   $Revision: 1.17.2.28 $  $Date: 2010/09/24 14:29:59 $
%   Copyright 2001-2010 The MathWorks, Inc.

% We will grow a number of cell strings as we generate the code. Hence we
% suppres this warning.
%#ok<*AGROW>

dsdb = getdsdb;
fitdb = getfitdb;
if isempty(down(dsdb)) && isempty(down(fitdb))
    errordlg( ...
        'Cannot generate file when no data sets or fits exist.', ...
        'Error Generating File', ...
        'modal' );
    return
end

% Parse inputs looking for a handle to the CFTOOL figure and where to write the
% results to
[cffig, fcnname, displayCodeFcn] = iParseInputs( varargin{:} );

% Set up some variables for later
allprop = {'Color' 'Marker' 'LineStyle' 'LineWidth' 'MarkerSize'};
dslist = cell(0,5);
ax = findall(cffig,'Type','axes','Tag','main');
xlim = get( ax, 'XLim' );
leginfo = {};
rleginfo = {};
showlegend = isequal(cfgetset('showlegend'),'on');
if showlegend
    % The user has selected to show the legend, still the possibility
    % exists that there is no legend. So we'll see if there is a legend
    % attached to the axis.
    legh = legend(ax);
    showlegend = ~isempty( legh );
end
if showlegend
    leginfo = cfgetlegendinfo(legh);
    mainlegloc = get(legh,'Location');
    uic1 = get(legh,'parent');          % uicontainer holding legend
    if isequal(mainlegloc,'none')
        oldu = get(legh,'units');
        legpos = get(legh,'position');     % legend position
        legpos = hgconvertunits(cffig,legpos,oldu,'normalized',uic1);
        mainlegloc = legpos(1:2);
    end
    ax = findall(cffig,'Type','axes','Tag','resid');
    residlegloc = 'none';
    if ~isempty(ax)
        legh = legend(ax);
        rleginfo = cfgetlegendinfo(legh);
        residlegloc = get(legh,'Location');
        uic2 = get(legh,'parent');          % uicontainer holding legend
        if isequal(residlegloc,'none')
            oldu = get(legh,'units');
            legpos = get(legh,'position');          % legend position
            legpos = hgconvertunits(cffig,legpos,oldu,'normalized',uic2);
            residlegloc = [1 .5] .* legpos(1:2);    % adjust to full figure
        end
        if isnumeric(mainlegloc)
            mainlegloc(2) = .5 + .5*mainlegloc(2); % adjust top legend too
        end
    end
end

% Create arrays to receive code text
blkc = cell(0,1);    % block of comment lines
blks = cell(0,1);    % block of setup lines
blkd = cell(0,1);    % block of data-related lines
blkf = cell(0,1);    % block of fit-related lines
blke = cell(0,1);    % block of lines at end

% Write introduction to data set section, including figure
% preparation code
blks{end+1} = iComment( xlate( 'Set up figure to receive data sets and fits' ) );
blks{end+1} = 'f_ = clf;';
blks{end+1} = 'figure(f_);';
fpos = cfgetfigurepos(cffig,'pixels');
blks{end+1} = sprintf('set(f_,''Units'',''Pixels'',''Position'',[%g %g %g %g]);',fpos);
if showlegend
    blks{end+1} = iComment( xlate( 'Line handles and text for the legend.' ) );
    blks{end+1} = 'legh_ = [];';
    blks{end+1} = 'legt_ = {};';
end
blks{end+1} = iComment( xlate( 'Limits of the x-axis.' ) );
blks{end+1} = 'xlim_ = [Inf -Inf];';

% Process each data set
nds = 0;
ds = down(dsdb);
arglist = {};
nonenames = 0;
while ~isempty(ds)
    nds = nds + 1;
    dsname = ds.name;
    dssrc = ds.source;
    
    % Was this data set created from raw data or from another source?
    smoothed = ~isempty(dssrc);
    if smoothed
        xname = dssrc{1,1};
        oldyname = dssrc{1,2};
        yname = sprintf('sm_.y%d',nds);
        wname = '';
    else
        xname = ds.xname;
        yname = ds.yname;
        wname = ds.weightname;
        oldyname = yname;
    end
    wtd = ~isempty(ds.weight) & ~isempty(wname);
    if isequal(lower(xname), cfGetNoneString() )
        xname = '';
    end
    
    % Each variable name becomes a function argument
    if ~isempty(xname) && isempty(strmatch(xname,arglist,'exact'))
        arglist{end+1} = xname; 
        newx = 1;
    else
        newx = 0;
    end
    if isempty(strmatch(oldyname,arglist,'exact')),
        arglist{end+1} = oldyname;
        newy = 1;
    else
        newy = 0;
    end
    if wtd && isempty(strmatch(wname,arglist,'exact'))
        arglist{end+1} = wname;
        neww = 1;
    else
        neww = 0;
    end
    
    % Create comment text associating data set with variable names
    blkc{end+1} =  ' ';
    blkc{end+1} = iComment( sprintf( 'Data from data set "%s":', dsname ) );
    if ~isempty(xname)
        blkc{end+1} = iComment( sprintf( '    X = %s:',xname ) );
    end
    blkc{end+1} =  iComment( sprintf( '    Y = %s:', oldyname ) );
    if wtd
        blkc{end+1} = iComment( sprintf( '    Weights = %s:', wname ) );
    else
        blkc{end+1} = iComment( xlate( '    Unweighted' ) );
    end
    
    % Create code to plot this data set into the figure we have created
    blkd{end+1} = ' ';
    blkd{end+1} = iComment( sprintf( '--- Plot data that was originally in data set "%s"',dsname) );
    if isempty(xname)
        nonenames = nonenames + 1;
        xname = sprintf('x_%d',nonenames);
        blkd{end+1} = sprintf('%s = (1:numel(%s))'';',xname,yname);
    elseif newx
        blkd{end+1} = sprintf('%s = %s(:);',xname,xname);
    end
    
    if newy,
        blkd{end+1} = sprintf('%s = %s(:);',oldyname,oldyname);
    end
    if neww, blkd{end+1} = sprintf('%s = %s(:);',wname,wname); end
    if smoothed
        % May need to reconstruct smoothed data from new input data
        for j=1:size(dssrc,1)
            span = dssrc{j,3};
            method = dssrc{j,4};
            degree = dssrc{j,5};
            blkd{end+1} = sprintf('%s = smooth(%s,%s,%g,''%s'',%g);',...
                yname,xname,oldyname,span,method,degree);
            oldyname = yname;
        end
    end
    dsline = ds.line;
    if ~isempty(dsline) && ishandle(dsline)
        propvals = get(dsline,allprop);
        [c,m,l,w,s] = deal(propvals{:});
        blkd{end+1} =  sprintf(...
            'h_ = line(%s,%s,''Parent'',ax_,''Color'',[%g %g %g],...', ...
            xname,yname,c(1),c(2),c(3));
        blkd{end+1} =  sprintf(...
            '     ''LineStyle'',''%s'', ''LineWidth'',%d,...', ...
            l,w);
        blkd{end+1} =  sprintf(...
            '     ''Marker'',''%s'', ''MarkerSize'',%d);', ...
            m,s);
        blkd{end+1} =  sprintf('xlim_(1) = min(xlim_(1),min(%s));',xname);
        blkd{end+1} =  sprintf('xlim_(2) = max(xlim_(2),max(%s));',xname);
        
        if showlegend
            blkd{end+1} = 'legh_(end+1) = h_;';
            blkd{end+1} = sprintf('legt_{end+1} = ''%s'';',quotedtext(dsname));
        end
    else
        blkd{end+1} = iComment( sprintf( ...
            'This data set does not appear on the plot. To add it to the plot, uncomment\nthe following lines and select the desired color and marker.' ) );
        propvals = '''Color'',''r'',''Marker'',''.'',''LineStyle'',''none''';
        blkd{end+1} = iComment( sprintf('   h_ = line(%s,%s,%s);', xname,yname,propvals) );
        blkd{end+1} = iComment( sprintf('   xlim_(1) = min(xlim_(1),min(%s));',xname) );
        blkd{end+1} = iComment( sprintf('   xlim_(2) = max(xlim_(2),max(%s));',xname) );
        
        if showlegend
            blkd{end+1} = iComment( '   legh_(end+1) = h_;' );
            blkd{end+1} = iComment( sprintf('   legt_{end+1} = ''%s'';',quotedtext(dsname)) );
        end
    end
    
    % Remember stuff about this data set
    dslist{nds,1} = dsname;
    dslist{nds,2} = xname;
    dslist{nds,3} = yname;
    if wtd, dslist{nds,4} = wname; end
    if ~isempty(dssrc), dslist{nds,5} = dssrc; end
    
    % Move to next data set
    ds = right(ds);
end

% Process each fit
nfit = 0;
ft = down(fitdb);
while ~isempty(ft)
    nfit = nfit + 1;
    ok = 1;
    if ~ft.isgood
        ok = 0;
    end
    if ok
        dsname = ft.dataset;
        fitname = ft.name;
        dsnum = strmatch(dsname,dslist(:,1),'exact');
        if isempty(dsnum)
            ok = 0;
        end
    end
    
    % altering indentation to avoid massive changes
    if ok
        xname = dslist{dsnum,2};
        yname = dslist{dsnum,3};
        wname = dslist{dsnum,4};
        wtd = ~isempty(wname);
        
        cf = ft.fit;
        ftype = type(cf);
        modeldef = ftype;
        
        % Create code to re-create this fit
        blkf{end+1} = ' '; % blank line
        blkf{end+1} = iComment( sprintf('--- Create fit "%s"',fitname) );
        fitargs = '';        % extra arguments for fit
        
        % If an excluded set exists in this fit, create an exclusion vector
        excluding = 0;
        if ~isequal(ft.outlier, cfGetNoneString() )
            outset = find(getoutlierdb,'name',ft.outlier);
            blkf{end+1} = ' '; % blank line
            blkf{end+1} = iComment( sprintf('Apply exclusion rule "%s"',ft.outlier) );
            
            % First deal with the exclusion vector
            evec = outset.exclude;
            n = length(evec);
            if n>0
                excluding = 1;
                
                % The generated code should throw an error if the exclusion rule
                % and the data are incompatible
                blkf{end+1} = sprintf('if length(%s)~=%d',xname,n);
                blkf{end+1} = iErrorCode( ...
                    'GenerateMFile:IncompatibleExclusionRule', ...
                    xlate( 'Exclusion rule ''%s'' is incompatible with ''%s''.' ), ...
                    ft.outlier, xname );
                blkf{end+1} = 'end';
                
                majority = (sum(evec) > .5*n);
                if majority
                    blkf{end+1} = sprintf('ex_ = true(length(%s),1);',xname);
                else
                    blkf{end+1} = sprintf('ex_ = false(length(%s),1);',xname);
                end
                if n>0
                    rowspec = makerowspec(evec == ~majority);
                    blkf{end+1} = sprintf('ex_(%s) = %d;',rowspec,~majority);
                end
            end
            
            % Apply restrictions based on X and Y
            if outset.restrictDomain
                comptxt1 = makecomptxt(xname,...
                    outset.domainLow,outset.domainLowLessEqual,...
                    outset.domainHigh,outset.domainHighLessEqual);
            else
                comptxt1 = '';
            end
            if outset.restrictRange
                comptxt2 = makecomptxt(yname,...
                    outset.rangeLow,outset.rangeLowLessEqual,...
                    outset.rangeHigh,outset.rangeHighLessEqual);
            else
                comptxt2 = '';
            end
            if isempty(comptxt1)
                comptxt1 = comptxt2;
            elseif ~isempty(comptxt2)
                comptxt1 = sprintf('%s | %s',comptxt1,comptxt2);
            end
            
            % Create a combined exclusion vector if necessary
            if ~isempty(comptxt1)
                excluding = 1;
                if n>0
                    blkf{end+1} = sprintf('ex_ = ex_ | %s;',comptxt1);
                else
                    blkf{end+1} = sprintf('ex_ = %s;',comptxt1);
                end
            end
        end
        
        % For custom models, get lots of information about the function
        ftypeargs = '';      % extra arguments for fittype
        start = [];          % starting points for iteration
        havefitopt = 0;
        if isempty(strmatch('custom',ftype))
            modeldef = sprintf('''%s''',modeldef);
        else
            % Get either the array of terms or the function expression
            if islinear(cf)
                modeldef = cell2text(linearterms(cf));
            else
                modeldef = cell2text(formula(cf));
            end
            
            % Get the variable and coefficient names
            cnames = coeffnames(cf);
            ftypeargs = sprintf(...
                ',...\n     ''dependent'',%s,''independent'',%s',...
                cell2text(dependnames(cf)), cell2text(indepnames(cf)));
            ftypeargs = sprintf('%s,...\n     ''coefficients'',%s',ftypeargs,...
                cell2text(cnames));
            
            % Get the problem parameter names, but this is not supported
            % in release one so it should be tested if they're added later
            pnames = probnames(cf);
            if ~isempty(pnames)
                ftypeargs = sprintf('%s,...\n     ''problem'',%s',ftypeargs,...
                    cell2text(pnames));
                pvals = num2cell(probvalues(cf));
                fitargs = sprintf('%s,''prob'',%s',fitargs,cell2text(pvals));
            end
        end
        
        % Get the fitoptions for this fit
        try
            oldopts = get(ft,'fitOptions');
        catch ignore %#ok<NASGU>
            oldopts = [];
        end
        
        % Use the same starting values as were used originally
        if (~isempty(oldopts)) && ~isempty(findprop(oldopts,'StartPoint'))
            start=oldopts.StartPoint;
        end
        
        % Examine the fitoptions object for this fit, and compare it
        % with a new one with default settings to figure out which
        % options to reproduce in the generated code
        ndiffs = 0;
        if isempty(oldopts)
            optfields = [];
        else
            oldmethod = oldopts.method;
            if isequal(oldmethod,'None') && isa(oldopts,'curvefit.smoothoptions')
                oldmethod = 'SmoothingSpline';
            end
            newopts = fitoptions('method',oldmethod);
            opttext = '';
            defaultfields = fields(newopts);
            optfields = fields(oldopts);
            
            % Remove fields that are done separately or should be ignored
            optfields(strcmpi('startpoint',optfields)) = [];
            optfields(strcmpi('method',optfields)) = [];
            optfields(strcmpi('exclude',optfields)) = [];
            optfields(strcmpi('display',optfields)) = [];
            optfields(strcmpi('weights',optfields)) = [];
        end
        for j=1:length(optfields)
            fname = optfields{j};
            k = find(strcmpi(fname,defaultfields));
            if length(k)==1
                optval = '';
                defval = '';
                try
                    optval = get(oldopts,fname);
                    if (strcmpi(fname,'lower') || strcmpi(fname,'upper')) ...
                            && all(isinf(optval))
                        optval = [];
                    end
                    defval = get(newopts,fname);
                catch ignore %#ok<NASGU>
                end
                if ~isequal(optval,defval) && (~isempty(optval) || ~isempty(defval))
                    ndiffs = ndiffs + 1;
                    opttext = sprintf('%s,''%s'',%s',opttext,fname,...
                        cell2text(optval));
                end
            end
        end
        if ndiffs>0
            blkf{end+1} = sprintf('fo_ = fitoptions(''method'',''%s''%s);',...
                oldmethod,opttext);
            fitargs = sprintf(',fo_%s',fitargs);
            havefitopt = 1;
        end
        
        if wtd
            blkf{end+1} =  sprintf('ok_ = isfinite(%s) & isfinite(%s) & isfinite(%s);',...
                xname,yname,wname);
        else
            blkf{end+1} =  sprintf('ok_ = isfinite(%s) & isfinite(%s);',xname,yname);
        end
        blkf{end+1} = 'if ~all( ok_ )';
        blkf{end+1} = iWarningCode( 'GenerateMFile:IgnoringNansAndInfs', ...
            xlate( 'Ignoring NaNs and Infs in data.' ) );
        blkf{end+1} = 'end';
        
        % Add things to either the fitoptions object or the fit call
        if ~isempty(start)
            blkf{end+1} = sprintf('st_ = [%s];',sprintf('%.20g ',start));
            if havefitopt
                blkf{end+1} = 'set(fo_,''Startpoint'',st_);';
            else
                fitargs = sprintf('%s,''Startpoint'',st_',fitargs);
            end
        end
        if wtd
            if havefitopt
                blkf{end+1} = sprintf('set(fo_,''Weight'',%s(ok_));',wname);
            else
                fitargs = sprintf('%s,''Weight'',%s(ok_)',fitargs,wname);
            end
        end
        if excluding
            if havefitopt
                blkf{end+1} = 'set(fo_,''Exclude'',ex_(ok_));';
            else
                fitargs = sprintf('%s,''Exclude'',ex_(ok_)',fitargs);
            end
        end
        
        blkf{end+1} = sprintf('ft_ = fittype(%s%s);',modeldef,ftypeargs);
        blkf{end+1} = ' '; % blank line
        blkf{end+1} = iComment( xlate( 'Fit this model using new data' ) );
        if excluding
            % The generated code should throw an error too many points are
            % excluded.
            blkf{end+1} = 'if sum(~ex_(ok_))<2';
            blkf{end+1} = iComment( xlate( 'Too many points excluded.' ) );
            blkf{end+1} = iErrorCode( ...
                'GenerateMFile:NotEnoughDataAfterExclusionRule', ...
                xlate( 'Not enough data left to fit ''%s'' after applying exclusion rule ''%s''.' ), ...
                fitname, ft.outlier );
            blkf{end+1} =  'else';
            blkf{end+1} =  sprintf('   cf_ = fit(%s(ok_),%s(ok_),ft_%s);',...
                xname,yname,fitargs);
            blkf{end+1} =  'end';
        else
            % If there is no exclusion rule, the generated code can just proceed
            % with the fit.
            blkf{end+1} =  sprintf('cf_ = fit(%s(ok_),%s(ok_),ft_%s);',...
                xname,yname,fitargs);
        end
        
        % For fits with simple coefficient representations, add a
        % helpful comment indicating what the estimates were
        categ = category(cf);
        cvals = [];
        try
            cvals = coeffvalues(cf);
        catch ignore %#ok<NASGU>
        end
        if ~isequal(categ,'spline') && ~isequal(categ,'interpolant') && ~isempty(cvals)
            blkf{end+1} = iComment( sprintf( ...
                'Alternatively uncomment the following lines to use coefficients from the\noriginal fit. You can use this choice to plot the original fit against new\ndata.' ) );
            blkf{end+1} = iComment( sprintf('   cv_ = %s;', cell2text(num2cell(cvals))) );
            blkf{end+1} = iComment( sprintf('   cf_ = cfit(ft_,cv_{:});') );
        end
        
        % Insert blank line before plotting code.
        blkf{end+1} = '';
        
        % Plot the fit and bounds if the original figure had them plotted
        if ~isempty(ft.line) && ishandle(ft.line)
            propvals = get(ft.line,allprop);
            [c,m,l,w,s] = deal(propvals{:});
            
            if iShowBoundsOn( ft )
                ptype = 'predobs';
                clev = ft.line.ConfLevel;
            else
                ptype = 'fit';
                clev = 0.95;
            end
            blkf{end+1} = iComment( xlate( 'Plot this fit' ) );
            blkf{end+1} = sprintf('h_ = plot(cf_,''%s'',%g);',ptype,clev);
            blkf{end+1} = sprintf( 'set(h_(1),''Color'',[%g %g %g],...', c(1),c(2),c(3));
            blkf{end+1} = sprintf('     ''LineStyle'',''%s'', ''LineWidth'',%d,...', l,w);
            blkf{end+1} = sprintf('     ''Marker'',''%s'', ''MarkerSize'',%d);', m,s);
            blkf{end+1} = iComment( xlate( 'Turn off legend created by plot method.' ) );
            blkf{end+1} = 'legend off;';
            if showlegend
                blkf{end+1} = iComment( xlate( 'Store line handle and fit name for legend.' ) );
                blkf{end+1} = 'legh_(end+1) = h_(1);';
                blkf{end+1} = sprintf('legt_{end+1} = ''%s'';', quotedtext(fitname));
            end
                        
            if isequal(ft.line.ShowBounds,'on')
                blkf{end+1} = 'if length(h_)>1';
                blkf{end+1} =  sprintf(...
                    '   set(h_(2:end),''Color'',[%g %g %g],...', ...
                    c(1),c(2),c(3));
                blkf{end+1} =  ...
                    '       ''LineStyle'','':'', ''LineWidth'',1,''Marker'',''none'');';
                if showlegend
                    blkf{end+1} = '   legh_(end+1) = h_(2);';
                    blkf{end+1} = sprintf('   legt_{end+1} = ''Pred bnds (%s)'';',...
                        quotedtext(fitname));
                end
                blkf{end+1} = 'end';
            end
            
            % Add residuals if the original figure had them
            if ~isempty(ft.rline) && ishandle(ft.rline)
                blkf{end+1} = ''; % blank line
                blkf{end+1} = iComment( xlate( 'Compute and plot residuals.' ) );
                propvals = get(ft.rline,allprop);
                [c,m,l,w,s] = deal(propvals{:});
                if excluding
                    blkf{end+1} =sprintf('res_ = %s(~ex_) - cf_(%s(~ex_));',...
                        yname,xname);
                    blkf{end+1} = sprintf('[x_,i_] = sort(%s(~ex_));',xname);
                else
                    blkf{end+1} = sprintf('res_ = %s - cf_(%s);',yname,xname);
                    blkf{end+1} = sprintf('[x_,i_] = sort(%s);',xname);
                end
                blkf{end+1} = 'axes(ax2_);';
                blkf{end+1} = 'hold on;';
                blkf{end+1} = sprintf(...
                    'h_ = line(x_,res_(i_),''Parent'',ax2_,''Color'',[%g %g %g],...',...
                    c(1),c(2),c(3));
                blkf{end+1} = sprintf(...
                    '     ''LineStyle'',''%s'', ''LineWidth'',%d,...', ...
                    l,w);
                blkf{end+1} = sprintf(...
                    '     ''Marker'',''%s'', ''MarkerSize'',%d);', ...
                    m,s);
                blkf{end+1} = 'axes(ax_);';
                blkf{end+1} = 'hold on;';
                if showlegend
                    blkf{end+1} = 'legrh_(end+1) = h_;';
                    blkf{end+1} = sprintf('legrt_{end+1} = ''%s'';',...
                        quotedtext(fitname));
                end
            end
        else
            blkf{end+1} = iComment( xlate( 'This fit does not appear on the plot' ) );
        end
    end
    % back to normal indentation
    
    % Move to next fit
    ft = right(ft);
end

% In setup section, create either one or two axes
% -- We setup two axes if we have residuals in CFTOOL
haveresiduals = ~strcmpi( cfgetset( 'residptype' ), 'none' );

if haveresiduals
    blks{end+1} = iComment( xlate( 'Axes for the main plot.' ) );
    blks{end+1} = 'ax_ = axes;';
    blks{end+1} = 'set(ax_,''Units'',''normalized'',''OuterPosition'',[0 .5 1 .5]);';
    blks{end+1} = iComment( xlate( 'Axes for the residuals plot.' ) );
    blks{end+1} = 'ax2_ = axes;';
    blks{end+1} = 'set(ax2_,''Units'',''normalized'',''OuterPosition'',[0 0 1 .5]);';
    blks{end+1} = 'set(ax2_,''Box'',''on'');';
    if showlegend
        blks{end+1} = iComment( xlate( 'Line handles and text for the residuals plot legend.' ) );
        blks{end+1} = 'legrh_ = [];';
        blks{end+1} = 'legrt_ = {};';
    end
else
    blks{end+1} = iComment( xlate( 'Axes for the plot.' ) );
    blks{end+1} = 'ax_ = axes;';
    blks{end+1} = 'set(ax_,''Units'',''normalized'',''OuterPosition'',[0 0 1 1]);';
end
blks{end+1} = 'set(ax_,''Box'',''on'');';
if isequal(cfgetset('showgrid'),'on')
    blks{end+1} = 'grid(ax_,''on'');';
    if haveresiduals
        blks{end+1} = 'grid(ax2_,''on'');';
    end
end
blks{end+1} = 'axes(ax_);';
blks{end+1} = 'hold on;';

% At end of data set section, set x axis limits
blkd{end+1} = ''; % blank line
blkd{end+1} = iComment( xlate( 'Nudge axis limits beyond data limits' ) );
blkd{end+1} = 'if all(isfinite(xlim_))';
blkd{end+1} = '   xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);';
blkd{end+1} = '   set(ax_,''XLim'',xlim_)';
if haveresiduals
    blkd{end+1} = '   set(ax2_,''XLim'',xlim_)';
end
% If the axis limits are infinite, then that means that we haven't
% plotted any data and the range will not be set up for the cfit/plot
% commands below.
blkd{end+1} = 'else';
blkd{end+1} = sprintf( '    set(ax_, ''XLim'',[%.20g, %.20g]);', xlim );
if haveresiduals
    blkd{end+1} = sprintf( '    set(ax2_,''XLim'',[%.20g, %.20g]);', xlim );
end
blkd{end+1} = 'end';


% Finish up
blke{end+1} = iComment( xlate( '--- Finished fitting and plotting data. Clean up.' ) ) ;
blke{end+1} = 'hold off;';
if showlegend
    blke{end+1} = iComment( xlate( 'Display legend' ) );
    blke{end+1} = sprintf('leginfo_ = %s; ', cell2text(leginfo));
    if isnumeric(mainlegloc)
		blke{end+1} = 'h_ = legend(ax_,legh_,legt_,leginfo_{:});';
        blke{end+1} = 'set(h_,''Units'',''normalized'');';
        blke{end+1} = 't_ = get(h_,''Position'');';
        blke{end+1} = sprintf('t_(1:2) = [%g,%g];',mainlegloc);
        blke{end+1} = 'set(h_,''Interpreter'',''none'',''Position'',t_);';
    else
        blke{end+1} = 'h_ = legend(ax_,legh_,legt_,leginfo_{:});';
        blke{end+1} = 'set(h_,''Interpreter'',''none'');';
    end
    
    if haveresiduals && nfit>0
        blke{end+1} = sprintf('leginfo_ = %s; % properties of resid legend', ...
            cell2text(rleginfo));
        if isnumeric(residlegloc)
            blke{end+1} = 'h_ = legend(ax2_,legrh_,legrt_,leginfo_{:});';
            blke{end+1} = 'set(h_,''Units'',''normalized'');';
            blke{end+1} = 't_ = get(h_,''Position'');';
            blke{end+1} = sprintf('t_(1:2) = [%g,%g];',residlegloc);
            blke{end+1} = 'set(h_,''Position'',t_);';
        else
            blke{end+1} = 'h_ = legend(ax2_,legrh_,legrt_,leginfo_{:});';
        end
        blke{end+1} = 'set(h_,''Interpreter'',''none'');';
    end
    blke{end+1} = iComment( xlate( 'Remove labels from x- and y-axes.' ) );
    blke{end+1} = 'xlabel(ax_,'''');';
    blke{end+1} = 'ylabel(ax_,'''');';
    if haveresiduals
        blke{end+1} = 'xlabel(ax2_,'''');';
        blke{end+1} = 'ylabel(ax2_,'''');';
        blke{end+1} = iComment( xlate( 'Add titles to the plots.' ) );
        blke{end+1} = 'title(ax_,''Data and Fits'');';
        blke{end+1} = 'title(ax2_,''Residuals'');';
    end
end

% Generate argument list text, i.e., the list of arguments to go into the
% "function" line of the file.
if isempty(arglist)
    argtext = '';
else
    argtext = sprintf('%s,',arglist{:});
    argtext = sprintf('(%s)',argtext(1:end-1));
end

% Make the help comment
blka = {};
blka{end+1} = sprintf( 'function %s%s', fcnname, argtext );
blka{end+1} = iH1Comment( fcnname, 'Create plot of data sets and fits' );
blka{end+1} = iHelpComment( sprintf( '%s%s', upper( fcnname ), upper( argtext ) ) );
blka{end+1} = iHelpComment( sprintf( 'Creates a plot, similar to the plot in the main Curve Fitting Tool,\nusing the data that you provide as input.  You can\nuse this function with the same data you used with CFTOOL\nor with different data.  You may want to edit the function to\ncustomize the code and this help message.' ) );
blka{end+1} = iHelpComment( '' ); % blank comment line
blka{end+1} = iHelpComment( sprintf( 'Number of data sets:  %d', nds ) );
blka{end+1} = iHelpComment( sprintf( 'Number of fits:  %d', nfit ) );

timestamp = iComment( sprintf( 'Auto-generated by MATLAB on %s', datestr( now ) ) );

% Put all of the code into a single string
mcode_str = '';
for j = 1:length( blka )
    mcode_str = sprintf( '%s%s\n', mcode_str, blka{j} );
end

for j=1:length(blkc)
    mcode_str = sprintf( '%s%s\n',mcode_str,blkc{j});
end

mcode_str = sprintf( '%s\n',mcode_str);
mcode_str = sprintf( '%s%s\n', mcode_str, timestamp );

mcode_str = sprintf( '%s\n',mcode_str);
for j=1:length(blks)
    mcode_str = sprintf( '%s%s\n',mcode_str,blks{j});
end

for j=1:length(blkd)
    mcode_str = sprintf( '%s%s\n',mcode_str,blkd{j});
end

for j=1:length(blkf)
    mcode_str = sprintf( '%s%s\n',mcode_str,blkf{j});
end

mcode_str = sprintf( '%s\n',mcode_str);
for j=1:length(blke)
    mcode_str = sprintf( '%s%s\n',mcode_str,blke{j});
end

displayCodeFcn( mcode_str );

% ------------------- display code in the MATLAB editor
function iDisplayCodeInEditor( mcode_str )

editorDoc = matlab.desktop.editor.newDocument(mcode_str);
editorDoc.smartIndentContents;
% Scroll document to line 1
editorDoc.goToLine(1);


% -------------------  write code to file
function iWriteCodeToFile(mcode_str,filename)

fid = fopen(filename,'w');
if(fid<0)
    emsg = sprintf( 'Could not create file: %s', filename );
    errordlg( emsg, 'Error Generating File', 'modal' );
    return
end
fprintf( fid, '%s',mcode_str);
fclose(fid);

% ------------------- double up quotes in text string
function a = quotedtext(b)
if ischar(b)
    a = strrep(b,'''','''''');
else
    a = sprintf('%.20g',b);
end

% ------------------- convert some text into a comment
function comment = iComment( text )
% This function is designed to handle multi-line comments as well as single line
% comments.
comment = iPrefix( '% ', text );

% ------------------- make an H1 comment
function comment = iH1Comment( functionName, description )
% An "H1 Comment" is the first line of the help comment. 
% 
% In code, an H1 comment is a single line consisting of
%   a comment symbol, followed by no space 
%   then the function name in upper case 
%   then a single space and 
%   then the (brief) description of the function.
comment = sprintf( '%%%s %s', upper( functionName ), description );

% ------------------- convert some text into a help comment
function comment = iHelpComment( text )
% A help comment has text indented three spaces from the comment symbol.
comment = iPrefix( '%   ', text );

% ------------------- add a prefix to some text
function text = iPrefix( prefix, text )
% Prefix some text with a given sting. This function is designed to handle
% multi-line text as well as single line text, i.e., the prefix is added to the
% start of each line of text.
%
% Example: A multi-line comment
%    iPrefix( '% ', sprintf( 'line one\nline two' );

NEWLINE = sprintf( '\n' );

% Add the prefix to the text
text = [prefix, text];

% Add the prefix after each new-line
text = strrep( text, NEWLINE, [NEWLINE, prefix] );

% ------------------- generate MATLAB code to display an error
function code = iErrorCode( id, msg, varargin )
% The variable arguments, varargin, must all be strings.
%
% The message, msg, and variable argument strings will have their contented
% "quoted" before being written out.
%
% See also QUOTEDTEXT.

% Ensure that the variable arguments are correctly quoted.
varargin = cellfun( @quotedtext, varargin, 'UniformOutput', false );

% Form a comma separated string from the variable arguments
otherArgs = sprintf( ', ''%s''', varargin{:} );
% Remove the initial ', '
otherArgs(1:2) = '';

% Generate the code for an error message -- spread across three lines: one line
% each for the id, the message and variable arguments
code = sprintf( 'error( ''%s'',...\n    ''%s'',...\n    %s );', ...
    id, quotedtext( msg ), otherArgs );

% ------------------- generate MATLAB code to display a warning
function code = iWarningCode( id, msg )
% The message, msg, string will have its contented "quoted" before being written
% out. 
%
% See also QUOTEDTEXT.

% Generate the code for a warning message -- spread across two lines: one line
% each for the id, and the message.
code = sprintf( 'warning( ''%s'',...\n    ''%s'' );', id, quotedtext( msg ) );

% ------------------- create text to re-create cell array
function a = cell2text( b )

% This is not a completely general-purpose routine, but it handles the sort
% of cell arrays used here.  A cell array containing a matrix, for
% instance, would not work
if ~iscell(b)
    if ischar(b)
        a = sprintf('''%s''',quotedtext(b));
    elseif length(b)==1
        a = sprintf('%.20g',b);
    else
        numtext = num2str(b,'%.20g ');
        if size(numtext,1)>1
            numtext = [numtext repmat(';',size(numtext,1),1)]';
            numtext = numtext(:)';
            numtext = numtext(1:end-1);
        end
        a = sprintf('[%s]',numtext);
    end
    return
end

if ~isempty(b)
    bj = b{1};
    if ischar(bj)
        a = sprintf('''%s''',quotedtext(bj));
    else
        a = sprintf(' %.20g',bj);
    end
    for j=2:length(b)
        bj = b{j};
        if ischar(bj)
            a = sprintf('%s, ''%s''',a,quotedtext(bj));
        elseif isscalar(bj)
            a = sprintf('%s, %.20g',a,bj);
        else
            a = sprintf('%s, [%s]',a,sprintf(' %.20g',bj));
        end
    end
else
    a = '';
end
a = sprintf('{%s}',a);

% ------------------- create text to generate a list of rows
function t = makerowspec(v)

% Find the start and end of each sequence
rownums = find(v(:));
t = '';

if ~isempty(rownums)
    diffs = diff(rownums);
    starts = rownums(logical([1;(diffs>1)]));
    ends   = rownums(logical([(diffs>1);1]));
    
    % Create some text to reproduce this sequence
    for j=1:length(starts)
        if starts(j) == ends(j)
            t = sprintf('%s %d',t,starts(j));
        elseif starts(j) == ends(j)-1
            t = sprintf('%s %d %d',t,starts(j),ends(j));
        else
            t = sprintf('%s (%d:%d)',t,starts(j),ends(j));
        end
    end
    t = t(2:end);
end

t = sprintf('[%s]',t);

% ------------------- create text to represent a comparison
function t = makecomptxt( varName, txtLow, isLowLessEqual, txtHigh, isHighLessEqual)

% for the isLowLessEqual and isHighLessEqual parameters,
%    --- "0" means "less than",
%    --- "1" means "less than or equal"

if isempty( txtLow ) || ~isfinite( str2double( txtLow ) )
    lowString = '';
elseif isLowLessEqual % == 1
    lowString = sprintf( '%s < %s', varName, txtLow );
else % isLowLessEqual == 0
    lowString = sprintf( '%s <= %s', varName, txtLow );
end

if isempty( txtHigh ) || ~isfinite( str2double( txtHigh ) )
    highString = '';
elseif isHighLessEqual % == 1
    highString = sprintf( '%s > %s', varName, txtHigh );
else % isHighLessEqual % == 0
    highString = sprintf( '%s >= %s', varName, txtHigh );
end

if isempty( lowString ) && isempty( highString )
    t = '';
elseif isempty( lowString )
    t = sprintf( '(%s)', highString );
elseif isempty( highString )
    t = sprintf( '(%s)', lowString );
else
    t = sprintf( '(%s | %s)', lowString, highString );
end

% ------------------- check for "show bounds" being on or off
function tf = iShowBoundsOn( ft )
% Display bounds is on if the line of the fit thinks that it is on and
% the type of fit is NOT a spline or an interpolant
ftype = category( ft.fit );
tf = isequal( ft.line.ShowBounds, 'on' ) ...
    && ~strcmpi( ftype, 'spline' ) ...
    && ~strcmpi( ftype, 'interpolant' );

% -------------------- Parse the inputs
function [cffig, fcnname, displayCodeFcn] = iParseInputs( cffig, filename )

if nargin < 1
    cffig = cfgetset('cffig'); 
end

% Have we been provided with the name of a file?
filenameGiven = nargin >= 2;

% We will wrtie to the editor if we have not been supplied with a file name
if ~filenameGiven
    % No file name supplied so make up a name of the function
    fcnname = 'createFit';
    % Request that output is sent to the editor
    displayCodeFcn = @iDisplayCodeInEditor;
else
    % Get file name with .m suffix, and get corresponding function name
    if length( filename ) < 2 || ~isequal( filename(end-1:end), '.m' )
        filename = sprintf( '%s.m', filename );
    end
    fcnname = filename(1:end-2);
    k = find( fcnname(1:end-1)=='\', 1, 'last' );
    if ~isempty( k )
        fcnname = fcnname(k+1:end);
    end
    k = find( fcnname(1:end-1)=='/', 1, 'last' );
    if ~isempty( k )
        fcnname = fcnname(k+1:end);
    end
    
    % Request that output is sent to the appropriate file
    displayCodeFcn = @(str) iWriteCodeToFile( str, filename );
end
