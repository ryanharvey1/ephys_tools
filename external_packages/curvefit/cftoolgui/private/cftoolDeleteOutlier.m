function cftoolDeleteOutlier(name)
%CFTOOLDELETEOUTLIER   Delete an exclusion rule from the outlier database
%
%   CFTOOLDELETEOUTLIER(NAME)

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2007/05/23 18:37:17 $ 

outdb = getoutlierdb;
h = find( outdb, 'name', name );
delete( h );
