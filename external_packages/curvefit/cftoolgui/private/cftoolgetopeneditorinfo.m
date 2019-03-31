function [name, isnormal, newfit, type, dataset, outlier, hint, results, usedefaultsmoothingparam] ...
                = cftoolgetopeneditorinfo(cftoolFit)
%CFTOOLGETOPENEDITORINFO Get information needed to open the cftool fit editor

%   $Revision: 1.1.6.3 $  $Date: 2007/06/14 04:55:09 $
%   Copyright 2004-2007 The MathWorks, Inc.

cftoolFit=handle(cftoolFit);

name = cftoolFit.name;
fitOptions = cftoolFit.fitOptions;
isnormal = false;
if isequal(fitOptions.Normalize, 'on')
    isnormal = true;
end
usedefaultsmoothingparam = false;
type = cftoolFit.type;
if isempty(type)
    newfit=true;
    dataset='';
    outlier='';
    hint='';
    results='';
else
    newfit = false;
    dataset = cftoolFit.dataset;
    outlier = cftoolFit.outlier;
    hint = cftoolFit.hint;
    results = cftoolFit.results;
    try
        if strcmp(type,'Smoothing Spline') && ...
                ~isempty(cftoolFit.output.p) && ...
                isempty(fitOptions.SmoothingParam)
            hint = sprintf('%g', cftoolFit.output.p);
            usedefaultsmoothingparam = true;
        end
    catch
       % Leave hint as originally assigned
    end
end
    


