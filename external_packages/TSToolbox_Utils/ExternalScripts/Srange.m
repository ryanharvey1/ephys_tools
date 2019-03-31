function y = Srange(x)
%RANGE  The range is the difference between the maximum and minimum values. 
%   Y = RANGE(X) calculates the range of the input.
%   For matrices RANGE(X) is a vector containing the range for each column.

%   Copyright 1993-2000 The MathWorks, Inc. 
%   $Revision: 2.7 $  $Date: 2000/05/26 18:53:32 $

y = max(x) - min(x);
