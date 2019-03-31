function values = previewValues(data)
% previewValues returns values for plotting
%
% previewValues(data) returns a (1,3) cell array. DATA is sftoolgui.Data.
% If there is a mismatch in sizes, the cell array has empty values. If all
% specified elements are equal, the cell array contains the actual values
% for those items that are specified and N 0s for those items not
% specified, where N is the number of elements

%   Copyright 2011 The MathWorks, Inc.
%   $Revision: 1.1.4.1 $    $Date: 2011/07/18 00:31:06 $

values = cell(1,3);
if areNumSpecifiedElementsEqual(data)
    if isCurveDataSpecified(data)
        [x, y] = getCurveValues(data);
        z = [];
    else
        [x, y, z] = getValues(data);
    end
    values{1} = x;
    values{2} = y;
    values{3} = z;
    
    isEmpty = cellfun( @(c) isempty( c ), values);
    n = cellfun( @(c) numel( c ), values(~isEmpty) );
    
    if ~isempty(n)
        for i = find( isEmpty )
            values{i} = zeros( 1, n(1) );
        end
    end
end
end

