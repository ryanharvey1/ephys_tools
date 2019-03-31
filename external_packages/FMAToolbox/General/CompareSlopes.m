%CompareSlopes - Perform linear regression on two groups and compare slopes.
%
%  USAGE
%
%    [p,t,df,slope1,slope2] = CompareSlopes(x1,y1,x2,y2)
%
%    x1,y1          data for group 1
%    x2,y2          data for group 2
%
%  OUTPUT
%
%    p              test probability (e.g. p<0.05 => significant)
%    t              t-test statistics
%    df             degrees of freedom
%    slope1,slope2  regression slopes

% Copyright (C) 2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [p,t,df,slope1,slope2] = CompareSlopes(x1,y1,x2,y2)

% Compute linear regressions
s1 = regstats(y1(:),x1(:),1,'tstat');
s2 = regstats(y2(:),x2(:),1,'tstat');

% Get N, slopes and SEM
n1 = length(x1);
n2 = length(x2);
slope1 = s1.tstat.beta;
slope2 = s2.tstat.beta;
sem1 = s1.tstat.se;
sem2 = s2.tstat.se;

% T-test
df = n1+n2-4;
t = (slope1-slope2)/sqrt(n1*sem1^2+n2*sem2^2);
p = 2*tcdf(-abs(t),df);






%  x = [x1(:);x2(:)];
%  y = [y1(:);y2(:)];
%
%  group = [zeros(length(x1),1);ones(length(x2),1)];
%
%  X = [group x x.*group];
%
%  s = regstats(y,X);
%
%  r =  [s.tstat.beta s.tstat.se s.tstat.t s.tstat.pval];




%  [slope1,u1,u2,u3,stats1] = regress(y1(:),x1(:));
%  [slope2,u1,u2,u3,stats2] = regress(y2(:),x2(:));
%
%  n1 = length(x1);
%  sse1 = stats1(end)/n1;
%  n2 = length(x2);
%  sse2 = stats2(end)/n2;
%
%  t = (slope1-slope2)/sqrt(sse1+sse2);
%  p = 2*tcdf(-abs(t),n1+n2-4);
%
%  r = [slope1 slope2 t p];



%  x = [x1(:);x2(:)];
%  y = [y1(:);y2(:)];
%
%  group = [zeros(length(x1),1);ones(length(x2),1)];
%
%  X = [group x x.*group];
%
%  [slope2,u1,u2,u3,stats2] = regress(y,X);
%  slope2, stats2
