function [V,W] = InterpolateTrackerZeros(X,Y)
%
% [V,W] = InterpolateTrackerZeros(X,Y)
% 
% Interpolate linearly the blackout regions (blocks) where the tracker 
% lost track of animal and returned x=y=0.   
%  
%  inputs:
%  X,Y are(c)tsd objects returned by [X,Y]=LoadPosition('VT1.pascii');
%
%  outputs:
%  V,W are corrected tsd objects; the missing position coordinates are 
%  linearly interpolated between last know position before a blackout block
%  and first known position after the blackout block
%
% PL 2000
% Version 0.2


x = Data(X);
y = Data(Y);

i0 = find(x+y == 0);
if length(i0) == 0
   % no corrections necessary
   V=X;
   W=Y;
   return;
end%if

ndata = length(x);

% check if last point is zero; if yes, then truncate all trailing zeros!
while(x(ndata) == 0 & y(ndata) == 0)
   x(ndata)=[];
   y(ndata)=[];
   ndata = length(x);
end%if
 
% interpolate linearly each block of consecutive zeros
for i = 1:ndata
   if(x(i) == 0 & y(i) == 0)
      nzeros = 1;
      ibstart = i-1;
      i=i+1;
      while (x(i) == 0 & y(i) == 0)   
         % count consecutive zeros
         nzeros = nzeros + 1;
         i=i+1;
      end%while
      ibend = i;
      x(ibstart:ibend)=linspace(x(ibstart),x(ibend),nzeros+2);
      y(ibstart:ibend)=linspace(y(ibstart),y(ibend),nzeros+2);
   end%if
end%for

ts = Range(X,'ts');
V = tsd(ts(1:ndata),x);
W = tsd(ts(1:ndata),y);


