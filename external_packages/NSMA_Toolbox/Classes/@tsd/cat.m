function tsdOUT = cat(varargin)

% tsd/cat  Concatenates list of n (c)tsd objects into a single tsd object
%
% tsdOUT = cat(tsd1, tsd2, ..., tsdn)
%
% INPUTS: 
%       tsd1, tsd2, ... tsdn - each one is either a ctsd or tsd 
% OUTPUTS:
%       tsdOUT - a tsd (not ctsd) that is a concatenation of all the inputs.
%
% ADR 1998, version L4.1, last modified '98 by ADR

% RELEASED as part of MClust 2.0
% See standard disclaimer in ../Contents.m

if ~isa(varargin{1}, 'tsd') & ~isa(varargin{1}, 'ctsd')
   error(['Initial tsd is not of type "[c]tsd".']);
end

tsdOUT = tsd(Range(varargin{1},'ts'), full(Data(varargin{1})));

for iTSD = 2:length(varargin)
   % First, all inputs must be tsd or ctsd 
   if ~isa(varargin{iTSD}, 'tsd') & ~isa(varargin{iTSD}, 'ctsd')
      error(['Input ', num2str(iTSD), 'is not a "[c]tsd"']);
   end
   % then all Data must be same dimension in non-time D
   szOUT = size(Data(tsdOUT));
   szTSD = size(Data(varargin{iTSD}));   
   if szOUT(2:length(szOUT)) ~= szTSD(2:length(szTSD))
      error(['Data size mismatch: input ', num2str(iTSD), '.']);
   end   
   % check to make sure times ok
   if StartTime(varargin{iTSD}) < EndTime(tsdOUT)
      error(['Time mismatch: input ', num2str(iTSD), 'starts before previous data ends.']);
   end

   tsdOUT.t = cat(1,tsdOUT.t, Range(varargin{iTSD},'ts'));
   tsdOUT.data = cat(1,tsdOUT.data, full(Data(varargin{iTSD})));

end