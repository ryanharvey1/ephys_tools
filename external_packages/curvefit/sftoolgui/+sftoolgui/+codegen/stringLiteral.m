function b = stringLiteral( a )
% STRINGLITERAL   Convert a string to a string literal for use in code
%
%   B = STRINGLITERAL( A ) converts the string A into the string
%   literal B. B will be valid for insertion in MATLAB code.

%   Copyright 2008-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2009/03/20 18:17:33 $

% Replace any single quote with double quotes
b = strrep( a, '''', '''''' );

end
