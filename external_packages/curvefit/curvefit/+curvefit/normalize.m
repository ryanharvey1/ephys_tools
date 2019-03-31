function [ normalizedData, meanData, stdData ] = normalize( data )
% NORMALIZE Normalize data
%
%   [NORMALIZEDDATA, MEANDATA, STDDATA] = NORMALIZE( DATA ) normalizes
%   DATA. NORMALIZEDDATA is the normalized data. MEANDATA is the mean of
%   DATA. STDDATA is the standard deviation of the data.

%   Copyright 2010 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2011/03/22 18:29:37 $ 

    meanData =  mean(data);
    
    stdData = std( data );
    % If the standard deviation is zero, then don't scale the data.
    stdData(stdData==0) = 1;
    
    % xdata = (xdata - meanx)/stdx with scalar expansion
    normalizedData = bsxfun( @rdivide, bsxfun( @minus, data, meanData ), stdData );
end

