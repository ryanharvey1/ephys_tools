function postConstructorSetup(h)
%POSTCONSTRUCTORSETUP   Set property values after construction.
%
%   POSTCONSTRUCTORSETUP(OBJ)

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2010/05/10 16:59:51 $ 

h.LowessOptionsVersion = 1;
h.method = 'LowessFit';
h.Robust = 'off';
h.Span = 0.25;
end
