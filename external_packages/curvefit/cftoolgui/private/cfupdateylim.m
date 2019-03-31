function cfupdateylim
%CFUPDATEYLIM Update the y axis min/max values

%   $Revision: 1.7.2.7 $  $Date: 2007/06/14 04:55:12 $
%   Copyright 2001-2007 The MathWorks, Inc.

dminmax = [];                    % to indicate y data limits

% Check any datasets with a plotting flag on
a = cfgetalldatasets;
fitdb = getfitdb;
for j=1:length(a)
   b = a{j};
   if b.plot == 1
      dminmax = combineminmax(dminmax,b.ylim);
   else
      flist = find(fitdb,'dataset',b.name,'plot',1);
      if ~isempty(flist)
         dminmax = combineminmax(dminmax,b.ylim);
      end
   end
end

% Check y limits of all fits
havedata = ~isempty(dminmax);
pminmax = dminmax;              % limits to use for plot
dy = diff(dminmax);
a = cfgetallfits;
for j=1:length(a)
   b = a{j};
   if b.plot==1 && ~isempty(b.line) && ishandle(b.line)
      fminmax = get(handle(b.line),'YLim');    %function limits
      if havedata && ~isempty(fminmax)
         if fminmax(2)<=dminmax(2)+dy && fminmax(1)>=dminmax(1)-dy
            % Expand range to include fit if extrapolation is reasonable
            pminmax(1) = min(pminmax(1),fminmax(1));
            pminmax(2) = max(pminmax(2),fminmax(2));
         else
            % Otherwise study the fit to see how much of it to plot
            ydata = sort(get(b.line,'ydata'));
            ydata(isnan(ydata)) = [];
            n = length(ydata);
            if n>4
               n4 = ceil(n/4);
               q1 = ydata(n4);
               q3 = ydata(n+1-n4);
               iqr = max(q3-q1,dy);
               fmin = max(ydata(1), q1-iqr);
               fmax = min(ydata(end), q3+iqr);
               pminmax(1) = min(pminmax(1),fmin);
               pminmax(2) = max(pminmax(2),fmax);
            end
         end
      else
         pminmax = combineminmax(pminmax,fminmax);
      end
   end
end

% Update main axis y limits
cfSetAxesLimits( 'YLim', pminmax );

% Update residual y axis if one exists
cffig = cfgetset('cffig');
ax2 = findall(cffig,'Type','axes','Tag','resid');
if ~isempty(ax2)
   set(ax2,'ylimmode','auto');
   ylimits = get(ax2,'ylim');
   set(ax2,'ylimmode','manual');
end

% ------------ Helper to combine old and new minmax values
function bothmm = combineminmax(oldmm,newmm)

if isempty(oldmm)
   bothmm = newmm;
elseif isempty(newmm)
   bothmm = oldmm;
else
   bothmm = [min(oldmm(1),newmm(1)) max(oldmm(2),newmm(2))];
end
