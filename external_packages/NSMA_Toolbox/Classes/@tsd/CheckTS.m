function bool = CheckTS(varargin)

% tsd/CheckTS  Checks whether all timestamps in list of (c)tsd objects are identical
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

% RELEASED as part of MClust 2.0
% See standard disclaimer in ../Contents.m


R0 = Range(varargin{1},'ts');
for iX = 2:length(varargin)
   R1 = Range(varargin{iX},'ts');
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