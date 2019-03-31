function [name, dL, dH, dLLE, dHLE, rL, rH, rLLE, rHLE, ds, ex, X, Y] = cftoolgetcopyexcludeinfo(outlier)
%CFTOOLGETCOPYEXCLUDEINFO Get information needed to copy a cftool exclusion rule 

%   $Revision: 1.1.6.3 $  $Date: 2007/11/09 19:48:50 $
%   Copyright 2004-2007 The MathWorks, Inc.

outlier = handle(outlier);

name = outlier.name;

dL = outlier.domainLow;
dH = outlier.domainHigh;
dLLE = outlier.domainLowLessEqual;
dHLE = outlier.domainHighLessEqual;

rL = outlier.rangeLow;
rH = outlier.rangeHigh;
rLLE = outlier.rangeLowLessEqual;
rHLE = outlier.rangeHighLessEqual;

ds = outlier.dataset;
ex = outlier.exclude;

dataset = i_datasetFromName( ds );
if isempty( dataset ),
    X = [];
    Y = [];
else
    [X, Y] = cfviewdata( dataset );
end

end

function dataset = i_datasetFromName( name )
if isempty( name )
    dataset = [];
else
    dsdb = getdsdb;
    dataset = dsdb.down;
    while ~isempty( dataset ) && ~isequal( dataset.name, name )
        dataset = dataset.right;
    end
end
end
