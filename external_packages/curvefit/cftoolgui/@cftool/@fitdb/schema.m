function schema
%
% $Revision: 1.8.2.2 $
% Copyright 2000-2004 The MathWorks, Inc.


pk = findpackage('cftool');

% Create a new class

%   $Revision: 1.8.2.2 $  $Date: 2005/03/07 17:24:59 $
c = schema.class(pk, 'fitdb');

schema.prop(c, 'current', 'string');

p=schema.prop(c, 'listeners', 'MATLAB array');
p.AccessFlags.Serialize = 'off';

p=schema.prop(c, 'newOptions', 'handle');
p.AccessFlags.Serialize = 'off';

p=schema.prop(c, 'newModel', 'string');
p.AccessFlags.Serialize = 'off';

p=schema.prop(c, 'newCoeff', 'string vector');
p.AccessFlags.Serialize = 'off';

p=schema.prop(c, 'newProps', 'string vector');
p.AccessFlags.Serialize = 'off';
