function schema
% A helper function for CFTOOL

% $Revision: 1.28.2.6 $  $Date: 2009/05/07 18:17:41 $
% Copyright 2000-2009 The MathWorks, Inc.

pk = findpackage('cftool');

% Create a new class called fit
c = schema.class(pk, 'fit');

% Add properties

% the fit itself (a cfit object)
p=schema.prop(c, 'fit', 'MATLAB array');
p.AccessFlags.Serialize = 'off';

% Equation Name -- This is a virtual property to get the
p = schema.prop( c, 'equationname', 'string' );
p.AccessFlags.Serialize  = 'off';
p.AccessFlags.PublicSet  = 'off';
p.AccessFlags.PrivateSet = 'off';
p.GetFunction = @iGetEquationName;

% Is Good -- did the fit complete without error
p=schema.prop(c, 'isGood', 'bool');

% the dataset for this fit
p=schema.prop(c, 'dshandle', 'MATLAB array');
p.AccessFlags.Serialize='off';

% the options for this fit
p=schema.prop(c, 'fitOptions', 'handle');

% labels
schema.prop(c, 'name', 'string');

schema.prop(c, 'type', 'string');
schema.prop(c, 'dataset', 'string');
schema.prop(c, 'results', 'string');
p=schema.prop(c, 'outlier', 'string');

% goodness
p=schema.prop(c, 'goodness', 'MATLAB array');
p=schema.prop(c, 'sse', 'double');
p=schema.prop(c, 'rsquare', 'double');
p=schema.prop(c, 'adjsquare', 'double');
p=schema.prop(c, 'rmse', 'double');
p=schema.prop(c, 'dfe', 'double');
p=schema.prop(c, 'ncoeff', 'double');

% output
p=schema.prop(c, 'output', 'MATLAB array');

% ingredients for confidence bounds
schema.prop(c, 'R', 'MATLAB array');

% a place to store plotting state
p=schema.prop(c, 'plot', 'int32');
p.AccessFlags.Serialize='off';
p=schema.prop(c, 'line', 'MATLAB array');
p.AccessFlags.Serialize='off';
p=schema.prop(c, 'rline', 'MATLAB array');   % residual line
p.AccessFlags.Serialize='off';

% A hint to make it easier to open in the editor.  
schema.prop(c, 'hint', 'string');

% starting values
schema.prop(c, 'start', 'MATLAB array');

% x and y limits
schema.prop(c, 'xlim', 'MATLAB array');
schema.prop(c, 'ylim', 'MATLAB array');

% a place to store listeners
p=schema.prop(c, 'listeners', 'MATLAB array');
p.AccessFlags.Serialize = 'off';

% Remember line color, marker, style, etc.
p=schema.prop(c, 'ColorMarkerLine', 'MATLAB array');

%-----------------------------------------------------------------------
function val = iGetEquationName( obj, ignore ) 
% "Get" function for equation name propterty
if isempty( obj.fit ),
    val = cfswitchyard( 'cfGetNoneString' );
else
    val = prettyname( obj.fit );
    if strcmpi( val, 'Model' ),
        % Let's assume it's a custom equation and use the hint as the name
        val = obj.hint;
    end
end
