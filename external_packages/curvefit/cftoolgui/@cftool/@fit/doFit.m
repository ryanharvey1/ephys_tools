function [results, cftoolFit] = doFit( cftoolFit, equation )
%DOFIT  Fit a fit-object to data
%
%   DOFIT(FIT, EQUATION)

%   Copyright 2001-2008 The MathWorks, Inc.
%   $Revision: 1.1.8.3 $    $Date: 2008/07/18 18:36:01 $


dataset=get(cftoolFit, 'dshandle');
fitOptions=get(cftoolFit, 'fitOptions');

% clear everything out
cftoolFit.plot=0;
convmsg = [];

% Get the equation name form the fit object
% convert the custom method string into the corresponding fittype object
if length(equation) > 7 && isequal(equation(1:8),'custom: ')
    equation = equation(9:end);
    equation = cfswitchyard( 'managecustom', 'get', equation );
    % Note that this returns a FITTYPE object
end

%if this method has a display, turn it off
if ~isempty(findprop(fitOptions,'Display'))
    fitOptions.Display='off';
end

%---------------------------------
% DO THE FIT
%---------------------------------
try
    [f, g, out, warnstr, errstr, convmsg] = fit( dataset.X, dataset.Y, ...
        equation, fitOptions );
catch e
    f = [];
    errstr= e.message;
    warnstr = '';
    out = [];
    g = [];
end

% Store all the information about this fit back into the cftool.fit object.
cftoolFit.fit=f;
if isempty( f ),
    % The fit failed in some way so we use the given fit options. We
    % make an explicit copy of them so that we don't end up with a
    % reference that some other object is also using.
    cftoolFit.fitOptions = copy( fitOptions );
else
    % We want to ensure that we get the fit options that were actually
    % used in the fit. Therefore we get them from the cfit object.
    % Another benefit of getting this way is that we get a proper copy
    % of the fit options and don't end up with multiple objects pointing
    % to the same fit options.
    cftoolFit.fitOptions = fitoptions( f );
end
cftoolFit.goodness=g;
cftoolFit.output=out;

% reset goodness of fit measures
cftoolFit.sse=NaN;
cftoolFit.rsquare=NaN;
cftoolFit.dfe=NaN;
cftoolFit.adjsquare=NaN;
cftoolFit.rmse=NaN;
cftoolFit.ncoeff=NaN;

% Check to see if a fit was created
if ~isempty(f)
    % A fit was created, this is good, :)
    cftoolFit.isGood=true;
    cftoolFit.sse=g.sse;
    cftoolFit.rsquare=g.rsquare;
    if isfield(g,'dfe')
        cftoolFit.dfe=g.dfe;
    else
        cftoolFit.dfe=0;
    end
    cftoolFit.adjsquare=g.adjrsquare;
    cftoolFit.rmse=g.rmse;
    if ~isempty(out.numparam)
        cftoolFit.ncoeff=out.numparam;
    end
    % Plot the line
    set(cftoolFit.line,'String',cftoolFit.name);
    cftoolFit.plot=1;

    % Generate the results text and store it in the fit
    clev = cfgetset( 'conflev' );
    resultsCellstr = genresults( f, g, out, warnstr, errstr, convmsg, clev );
    results = sprintf( '%s\n', resultsCellstr{:} );
    cftoolFit.results = results;
else
    % A fit was not created, this is not good, :(
    cftoolFit.isGood=false;
    
    % Because the fit failed, we don't have a CFIT object to store.
    % However, a FITTYPE is almost as good.
    cftoolFit.fit = fittype( equation );

    % The fit failed. However we can still generate the results string by using
    % any CFIT in place of the fit.
    resultsCellstr = genresults( cfit(), g, out, warnstr, errstr, convmsg );
    results = sprintf( '%s\n', resultsCellstr{:} );
    cftoolFit.results = results;

end

% CreateAFit uses hint to set up the state of the parameter panel when opening
% an existing fit. This hint is mostly set in CreateAFit, but we want to
% return the value of the smoothing parameter that was used in the fit,
% which may have been estimated by the fit.
if isequal( equation, 'smoothing' ) ...
        && isempty( fitOptions.smoothingparam ) ... % This means the default value was used
        && isstruct( out )
    % The fit command used the default smoothing parameter.  Stick
    % this back into the cftool.fit object's hint for later use.
    cftoolFit.hint = sprintf( '%g', out.p );
end
