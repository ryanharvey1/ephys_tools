function schema
%
% $Revision: 1.16.2.3 $	$Date: 2010/03/15 22:31:53 $
% Copyright 2000-2010 The MathWorks, Inc.


pk = findpackage('cftool');

% Create a new class called data

c = schema.class(pk, 'dataset');

% Add properties
p = schema.prop(c, 'x', 'MATLAB array');
p.GetFunction = @iEnsureFull;

p = schema.prop(c, 'y', 'MATLAB array');
p.GetFunction = @iEnsureFull;

p = schema.prop(c, 'weight', 'MATLAB array');
p.GetFunction = @iEnsureFull;

schema.prop(c, 'name', 'string');

schema.prop(c, 'xlength', 'double');

schema.prop(c, 'xname', 'string');
schema.prop(c, 'yname', 'string');
schema.prop(c, 'weightname', 'string');

schema.prop(c, 'xlim', 'MATLAB array');
schema.prop(c, 'ylim', 'MATLAB array');

p=schema.prop(c, 'plot', 'int32');
p.AccessFlags.Serialize = 'off';
p=schema.prop(c, 'line', 'MATLAB array');
p.AccessFlags.Serialize = 'off';

p=schema.prop(c, 'listeners', 'MATLAB array');
p.AccessFlags.Serialize = 'off';

schema.prop(c, 'ColorMarkerLine', 'MATLAB array');
schema.prop(c, 'source', 'MATLAB array');


function val = iEnsureFull(~, val)
% iEnsureFull( obj, val ) -- Get Function to ensure that the returned parameter
% value is full (not sparse)
val = full( val );
