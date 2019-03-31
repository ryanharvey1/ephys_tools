function [StartTS, EndTS] = SMatrix_StartEndTimes(S)
% return arrays of  start and end timesamps of a given Spike Matrix (S-Matrix)
%
% [StartTS, EndTS] = SMatrix_StartEndTimes(S)
%
% loops over all ts-objects in a given S-matrix and finds: 
%   SessStartTS  ... array of first timestamp of all cells (ts-objects)
%   SessEndTS    ... array of the last timestamps in all ts-objects
%
% PL Feb. 2003

% find SessStartTS, SessStopTS:
StartTS = zeros(size(S));
EndTS = zeros(size(S));
for i = 1:length(S)
    if length(data(S{i})) > 0             % count only cells with at least one spike
        StartTS(i) = StartTime(S{i});
        EndTS(i) = EndTime(S{i});  
    end
end
