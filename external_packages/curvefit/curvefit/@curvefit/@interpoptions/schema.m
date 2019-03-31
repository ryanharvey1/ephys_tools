function schema
% Schema for interp object.

% Copyright 2001-2005 The MathWorks, Inc.
% $Revision: 1.4.2.2 $  $Date: 2005/03/07 17:25:53 $

pk = findpackage('curvefit');

% Create a new class called interpoptions

schema.class(pk, 'interpoptions', pk.findclass('basefitoptions'));



