function [maxtab, mintab]=peakdetz(v, delta, lookformax, backwards)
%PEAKDET Detect peaks in a vector
%        [MAXTAB, MINTAB] = PEAKDETZ(V, DELTA, lookformax, backwards) finds 
%        the local maxima and minima ("peaks") in the vector V.
%        A point is considered a maximum peak if it has the maximal
%        value, and was preceded (to the left) by a value lower by
%        DELTA. MAXTAB and MINTAB consists of two columns. Column 1
%        contains indices in V, and column 2 the found values.
%
% Eli Billauer, 3.4.05 (Explicitly not copyrighted).
% This function is released to the public domain; Any use is allowed.
%
% ZN edit 04/2010: added option to specify looking for troughs or peaks
% first (lookformax variable: if 1, will look for peaks first, if 0 will
% look for troughs; default is look for peaks); and option to go backwards
% (so that find last instance of a peak/trough value instead of the first
% instance: backwards variable: if 1 will go backwards, if 0 or absent, 
% will go forwards); and changed it so that last min/max value will be 
% assigned

if nargin<3
    lookformax = 1;
end

maxtab = [];
mintab = [];

v = v(:); % Just in case this wasn't a proper vector

if (length(delta(:)))>1
  error('Input argument DELTA must be a scalar');
end

if delta <= 0
  error('Input argument DELTA must be positive');
end

if nargin<4 || backwards==0
    inc=1;
    first=1;
    last=length(v);
elseif backwards
    inc=-1;
    first=length(v);
    last=1;
end


mn = Inf; mx = -Inf;
mnpos = NaN; mxpos = NaN;

for ii=first:inc:last
  this = v(ii);
  if this > mx, mx = this; mxpos = ii; end
  if this < mn, mn = this; mnpos = ii; end
  
  if lookformax
    if this < mx-delta || (ii==last && ~isempty(mintab) && mx-delta>mintab(end,2))
      maxtab = [maxtab ; mxpos mx];
      mn = this; mnpos = ii;
      lookformax = 0;
    end  
  else
    if this > mn+delta || (ii==last && ~isempty(maxtab) && mn+delta<maxtab(end,2))
      mintab = [mintab ; mnpos mn];
      mx = this; mxpos = ii;
      lookformax = 1;
    end
  end
end

if isempty([mintab;maxtab])
    if lookformax
        if mx-mn>delta
            maxtab=[mxpos mx];
        end
    else
        if mx-mn>delta
            mintab=[mnpos mn];
        end
    end
end
