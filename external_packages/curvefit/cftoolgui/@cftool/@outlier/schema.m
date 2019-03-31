function schema
%
% $Revision: 1.6.2.2 $	$Date: 2004/12/06 16:33:33 $
% Copyright 1999-2004 The MathWorks, Inc.

pk = findpackage('cftool');

% Create a new class called outlier

c = schema.class(pk, 'outlier');

% Add properties
schema.prop(c, 'name', 'string');
p=schema.prop(c, 'dataset', 'string');

p=schema.prop(c, 'restrictDomain', 'bool');
p=schema.prop(c, 'domainLow', 'string');
p=schema.prop(c, 'domainHigh', 'string');
p=schema.prop(c, 'restrictRange', 'bool');
p=schema.prop(c, 'rangeLow', 'string');
p=schema.prop(c, 'rangeHigh', 'string');

% for these "equal" properties, "0" means "less than", "1" means "less than or equal"
p=schema.prop(c, 'domainLowLessEqual', 'double');
p=schema.prop(c, 'domainHighLessEqual', 'double');
p=schema.prop(c, 'rangeLowLessEqual', 'double');
p=schema.prop(c, 'rangeHighLessEqual', 'double');

p = schema.prop(c, 'exclude', 'MATLAB array');
p=schema.prop(c, 'length', 'double');

p=schema.prop(c, 'listeners', 'MATLAB array'); % place to store listeners
p.AccessFlags.Serialize = 'off';
