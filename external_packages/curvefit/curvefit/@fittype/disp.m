function disp(obj)
%DISP   DISP for FITTYPE.

%   Copyright 1999-2005 The MathWorks, Inc.
%   $Revision: 1.4.2.3 $  $Date: 2008/10/31 05:56:35 $



objectname = inputname(1);

[ignore,line2] = makedisplay(obj,objectname);

fprintf('     %s\n', line2);
