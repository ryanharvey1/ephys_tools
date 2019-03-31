function [S,E] = IntervalsAboveThreshold(v, threshold, MinIntervalLength , MinGapLength, plotWanted)
% 
% if v is a vector of a time series 
% find intervals (defined by start index S and end index E)  where v is 
% above threshold. Keep only intervals that are longer-equal than MinInteralLength and have gaps greater-equal MinGapLength
% (either one or both can be zero)
%
% if plotWanted ~= 0, a figure with v and the intervals is produced
% 

if length(v) < 4
    error(' Input vector must have at least 4 elements!');
end

[nr,nc] = size(v);
if(nc ~= 1)
    v = v';      % make sure v is a column vector
end

inth = find(v > threshold);
ii   = find(diff(inth) > 1);
start_i = [inth(1); inth(ii+1)];
end_i = [inth(ii); inth(end)];

S = start_i;
E = end_i;

% cut out gaps shorter than MinGapLength
i = 2;
while i < length(S)
  while (S(i) - E(i-1)) < MinGapLength & i < length(S)
    S(i) = [];
    E(i-1) = [];
  end
  i = i+1;
end

% cut out intervals shorter than MinIntervalLength
i=1;
while i <= length(S)
  while (E(i)-S(i)) < MinIntervalLength & i < length(S)
    S(i) = [];
    E(i) = [];
  end
  i = i+1;
end


% make figure for visual inspection
if plotWanted ~= 0
    figure;
    plot(v);
    axis tight;
    hold on;
    startend = [start_i end_i];
    plot(startend', threshold*ones(size(startend))');
    hold off;
end

