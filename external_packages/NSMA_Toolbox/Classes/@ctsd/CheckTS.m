function bool = CheckTS(varargin)

% ctsd/CheckTS  Checks whether all timestamps in list of (c)tsd objects are identical
%
% bool = CheckTS(X0, X1, X2, ...)
%
% INPUTS: 
%       X0, X1, X2... - each one is either a ctsd or tsd 
% OUTPUTS:
%       bool - 1 if timestamps all identical, 0 if not
%
% Works with combinations of ctsd and tsd
%
% ADR 1998, version L4.0, last modified '98 by ADR

% status: PROMOTED


% for speed, first check to make sure if they're all ctsd's...
flagAllCTSD = 1;
for iX = 1:length(varargin)
   if ~isa(varargin{iX}, 'ctsd')
      flagAllCTSD = 0;
   end
end

if flagAllCTSD
   % if they're all ctsd's, just check DT's and start/end times
   DT0 = DT(varargin{1});
   T0 = StartTime(varargin{1});
   T1 = EndTime(varargin{1});
   for iX = 2:length(varargin)
      if (DT0 ~= DT(varargin{iX}))
         bool = 0; 
         return
      elseif (T0 ~= StartTime(varargin{iX}))
         bool = 0;
         return
      elseif (T1 ~= EndTime(varargin{iX}))
         bool = 0;
         return
      end % if
   end % for
   bool = 1;
   return
else
   R0 = Range(varargin{1});
   for iX = 2:length(varargin)
      R1 = Range(varargin{iX});
      if (length(R0) ~= length(R1))
         bool = 0;                  % if not same length, can't be equal.
         return
      end
      if (min(R0 == R1) == 0)
         bool = 0;                  % if there are any non-equal elts, not equal
         return
      end
   end
   bool = 1;                        % nothing failed, must be ok
   return
end
