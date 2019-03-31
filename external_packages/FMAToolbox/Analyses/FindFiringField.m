%FindFiringField - Find firing field in a firing map.
%
% Find firing field in a firing map, i.e. the connex area around the peak,
% where the firing rates are above a given threshold.
%
% Note: This function should not be called directly. Rather, the firing field
% should be obtained with <a href="matlab:help FiringMap">FiringMap</a>, which in turn repeatedly calls
% FindFiringField with the appropriate parameters (repetitions are used to
% discard regions of spurious elevated firing rates).
%
%  USAGE
%
%    field = FindFiringField(x,y,M,N,threshold,map)
%
%    x              abscissa (column) of the peak
%    y              ordinate (row) of the peak
%    M              abscissae run from 1 to M
%    N              ordinates run from 1 to N
%    threshold      min firing rate within the field
%    map            firing map
%
%  SEE
%
%    See also DefineField, FiringMap, IsInField, PlotColorMap.

% Copyright (C) 2004-2006 by Michaël Zugaro
