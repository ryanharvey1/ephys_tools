function h = outlier(name, dataset, rd, dl, dh, rr, rl, rh, dle, dhe, rle, rhe, exclude, len)

% $Revision: 1.6.2.3 $  $Date: 2005/03/07 17:25:00 $
% Copyright 1999-2005 The MathWorks, Inc.

h = cftool.outlier;

if nargin==0
   h.dataset='';
   h.exclude=[];
   h.restrictDomain=false;
   h.domainLow='';
   h.domainHigh='';
   h.restrictRange=false;
   h.rangeLow='';
   h.rangeHigh='';
   h.domainLowLessEqual=0;
   h.domainHighLessEqual=0;
   h.rangeLowLessEqual=0;
   h.rangeHighLessEqual=0;
   h.length=0;
   h.name = '';
else
   h.dataset=dataset;
   
   if (len == 0)
      h.exclude=[];
   else
      h.exclude=exclude;
   end
   
   if (rd == 0)
      h.restrictDomain=false;
      h.domainLow='';
      h.domainHigh='';
   else
      h.restrictDomain=true;
      h.domainLow=dl;
      h.domainHigh=dh;
   end
   if (rr == 0)
      h.restrictRange=false;
      h.rangeLow='';
      h.rangeHigh='';
   else
      h.restrictRange=true;
      h.rangeLow=rl;
      h.rangeHigh=rh;
   end
   h.domainLowLessEqual=dle;
   h.domainHighLessEqual=dhe;
   h.rangeLowLessEqual=rle;
   h.rangeHighLessEqual=rhe;
   h.length=len;
   
   % assumes name is unique
   h.name = name;
end

% add it to the list of outliers
connect(h,getoutlierdb,'up');
cfgetset('dirty',true);   % session has changed since last save


list(2) = handle.listener(h,findprop(h,'name'),'PropertyPostSet',...
                          {@updatename});
list(1) = handle.listener(h,'ObjectBeingDestroyed', {@cleanup});
h.listeners=list;

%=============================================================================
function updatename(ignore1,ignore2)

cfgetset('dirty',true);   % session has changed since last save


%=============================================================================
function cleanup(ignore1,ignore2)

cfgetset('dirty',true);   % session has changed since last save

