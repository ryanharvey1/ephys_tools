function disp(obj)
%DISP   DISP for CFIT object.

%   Copyright 1999-2005 The MathWorks, Inc.
%   $Revision: 1.4.2.3 $  $Date: 2008/10/31 05:56:18 $

objectname = inputname(1);

[ignore,line2,line3,line4] = makedisplay(obj,objectname);

fprintf('     %s\n', line2);
fprintf('     %s\n', line3);
if ~isempty(line4)
    fprintf('     %s\n', line4); 
end
