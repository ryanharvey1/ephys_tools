function cfsmoothds(dsname,dataset,method,span,degree)
%CFSMOOTHDS

%   $Revision: 1.7.2.2 $  $Date: 2007/10/15 22:40:59 $
%   Copyright 2001-2007 The MathWorks, Inc.

%CFSMOOTHDS(NAME,DATASET,METHOD,SPAN,DEGREE)

ds = handle(dataset);
z = smooth(ds.x,ds.y,span,method,degree);

% Store information about the source of the data in this dataset
newsrc = {ds.xname ds.yname span method degree};
if ~isempty(ds.source)
   src = ds.source;
   src(end+1,:) = newsrc;
   newsrc = src;
end

cftool.dataset(ds.xname, ds.yname, cfGetNoneString() ,dsname,ds.x,z,[],newsrc);
