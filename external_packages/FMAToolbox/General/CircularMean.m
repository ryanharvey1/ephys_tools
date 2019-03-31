%CircularMean - Estimate the circular mean.
%
%  USAGE
%
%    [m,r] = CircularMean(angles,dim,centered)
%
%    angles         angles in radians
%    dim            optional dimension along which the mean should be computed
%    m              mean angle
%    r              resultant length
%
%  SEE
%
%    See also CircularVariance, CircularConfidenceIntervals, Concentration,
%    ConcentrationTest.

% Copyright (C) 2004-2008 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [m,r] = CircularMean(angles,dim)

if isempty(angles), m = []; return; end

if nargin < 2,
	dim = 1;
end

% isradians(angles);

angles = exp(i*angles);
try r = abs(nanmean(angles,dim)); catch;r=abs(mean(angles,dim));end
try m = angle(nanmean(angles,dim)); catch;m = angle(mean(angles,dim));end
