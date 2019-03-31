function tsa = tsd(t, Data)

% tsd  Creates tsd object (tsd = class of timestamped arrays of data)
%
% tsa = tsd(t, Data)
%
% INPUTS: 
%       t - list of timestamps - must be sequential, but don't have to be continuous
%       Data - data, possibly an array. If data is n-dimensional, then time should be the FIRST axis.
% OUTPUTS: 
%       tsa - a tsd object
%
% Completely compatible with ctsd.
%
% Methods
%    tsd/Range     - Timestamps used
%    tsd/Data      - Returns the data component
%    tsd/DT        - Returns the DT value (mean diff(timestamps))
%    tsd/StartTime - First timestamp
%    tsd/EndTime   - Last timestamp
%    tsd/Restrict  - Keep data within a certain range
%    tsd/CheckTS   - Makes sure that a set of tsd & ctsd objects have identical start and end times
%    tsd/cat       - Concatenate ctsd and tsd objects
%    tsd/Mask      - Make all non-mask values NaN
%
% ADR 1998, version L4.0, last modified '98 by ADR
% RELEASED as part of MClust 2.0
% See standard disclaimer in ../Contents.m

if nargin == 0
 tsa.t = NaN;
 tsa.data = NaN;
 tsa = class(tsa, 'tsd');
 return
end 

if nargin < 2
  error('tsd constructor must be called with T, Data');
end

tsa.t = t;
tsa.data = Data;
tsa = class(tsa, 'tsd');
