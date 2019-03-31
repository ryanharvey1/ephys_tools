function x = nanFromComplexElements( x )
% nanFromComplexElements   Replace each complex element, i.e., element with
% non-zero imaginary part, and set it to NaN.

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2011/07/18 00:30:42 $

if ~isreal( x )
    % Determine which elements of the array have non-zero imaginary part.
    tf = imag( x ) ~= 0;
    % Set those elements to NaN
    x(tf) = NaN;
end
end
