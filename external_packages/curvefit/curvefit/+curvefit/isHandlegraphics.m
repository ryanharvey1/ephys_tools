function tf = isHandlegraphics(obj)
%isHandlegraphics True for handlegraphics objects
%
%   tf = isHandlegraphics(OBJ) returns true if OBJ is a handlegraphics
%   object
%
%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date:  

tf = graphicsversion(obj, 'handlegraphics');
end
