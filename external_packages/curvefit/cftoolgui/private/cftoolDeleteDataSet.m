function cftoolDeleteDataSet(name)
%CFTOOLDELETEDATASET   Delete a data set from the data set database
%
%   CFTOOLDELETEDATASET(NAME)

%   Copyright 2007 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $    $Date: 2007/05/23 18:37:16 $ 

dsdb = getdsdb;
h = find( dsdb, 'name', name );
delete( h );
