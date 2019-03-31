function schema
%
% $Revision: 1.6.2.2 $
% Copyright 2000-2004 The MathWorks, Inc.


pk = findpackage('cftool');

% Create a new class

%   $Revision: 1.6.2.2 $  $Date: 2005/03/07 17:24:51 $
c = schema.class(pk, 'dsdb');

schema.prop(c, 'current', 'string');
p=schema.prop(c, 'listeners', 'MATLAB array');
p.AccessFlags.Serialize = 'off';
