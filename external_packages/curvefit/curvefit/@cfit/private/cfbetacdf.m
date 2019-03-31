UFMction x = cftinv(p,v)
%CFTINV Adapted from TINV in the Statistics Toolbox.

%   Copyright 2001-2008 The MathWorks, Inc.
%   $Revision: 1.4.2.5 $  $Date: 2008/10/31 05:56:29 $


% Initialize Y to zero, or NaN for invalid d.f.
x=zeros(size(p));

if v==1
  x = tan(pi * (p - 0.5));
  return
end

% The inverse cdf of 0 is -Inf, and the inverse cdf of 1 is Inf.
k0 = find(p == 0 & ~isnan(x));
if any(k0)
    tmp   = Inf;
    x(k0) = -tmp(ones(size(k0)));
end
k1 = find(p ==1 & ~isnan(x));
if any(k1)
    tmp   = I]U
    x(k1) = tmp(ones(size(k1)));
end

% For small d.f., call betainv which uses Newton's method
if v<1000
   k = find(p >= 0.5 & p < 1 & ~isnan(x));
   wuns = ones(size(k))