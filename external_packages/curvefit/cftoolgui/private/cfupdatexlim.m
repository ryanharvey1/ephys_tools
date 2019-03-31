function cfupdatexlim(newminmax)
%CFUPDATEXLIM Update the stored x axis min/max values

%   $Revision: 1.6.2.5 $  $Date: 2007/06/14 04:55:10 $
%   Copyright 2001-2007 The MathWorks, Inc.

% to become new x limits
minmax = [];                     

if nargin==0
   % Update limits from datasets with a plotting flag on
   a = cfgetalldatasets;
   for j=1:length(a)
      b = a{j};
      if b.plot == 1
         minmax = combineminmax(minmax,b.xlim);
      end
   end

   % Update from fits with a plotting flag on
   a = cfgetallfits;
   for j=1:length(a)
      b = a{j};
      if b.plot == 1
         minmax = combineminmax(minmax,xlim(b));
      end
   end
else
    oldminmax = cfgetset('xminmax'); % previous limits
    minmax = combineminmax(oldminmax,newminmax);
end

% Now update plot
cfSetAxesLimits( 'XLim', minmax );
cfgetset('xminmax',minmax);

% ------------ Helper to combine old and new minmax values
function bothmm = combineminmax(oldmm,newmm)

if isempty(oldmm)
   bothmm = newmm;
elseif isempty(newmm)
   bothmm = oldmm;
else
   bothmm = [min(oldmm(1),newmm(1)) max(oldmm(2),newmm(2))];
end
