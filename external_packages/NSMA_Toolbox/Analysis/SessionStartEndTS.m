function [SessStartTS, SessEndTS] = SessionStartEndTS(S)
% find start and end timesamps of a given Spike Matrix (S-Matrix)
%
% [SessStartTS, SessEndTS] = SessionStartEndTS(S)
%
% loops over all ts-objects in a given S-matrix and finds: 
%   SessStartTS  ... the approximate beginning of the session
%   SessEndTS    ... the approximate end of the session
% 
% To remove outliers, the mean and std of the start and end timestamps of each cell in S are computed
% and the following definitions are applied:
%   cutStartTS = mean(starttimes) - 3*std(starttimes);
%   cutEndTS   = mean(endtimes)   + 3*std(starttimes);
% Any cell with start or end timesamp outside the [cutStartTS, cutEndTS] interval is considered an outlier
% (outside the 99th percentile) and removed from the list. From the resulting list, the 
%   SessStartTS = min(startTS) and the 
%   SessEndTS = max(endTS)     are comupted.
%
%  This definition prevents an epoch to be dominated by spurious first or last timestamps of single t-files
% (which may be off by orders of magnitudes - which happens occasionally for unknown reasons).
%
% PL Feb. 2003

% find SessStartTS, SessStopTS:

[StartTS, EndTS] = SMatrix_StartEndTimes(S);

warning('off','MATLAB:divideByZero');
avgStartTS = nanmean(StartTS);
stdStartTS = nanstd(StartTS);
cutStartTS = avgStartTS - 3*stdStartTS;

avgEndTS = nanmean(EndTS);
stdEndTS = nanstd(EndTS);
cutEndTS = avgEndTS + 3*stdEndTS;
warning('off','MATLAB:divideByZero');

ixcut = [];
for i=1:length(StartTS)
    if StartTS(i) < cutStartTS
        ixcut(end+1) = i;
    end    
end
StartTS(ixcut) = [];

ixcut = [];
for i=1:length(EndTS)
    if EndTS(i) > cutEndTS
        ixcut(end+1) = i;
    end    
end
EndTS(ixcut) = [];

SessStartTS = min(StartTS);
SessEndTS = max(EndTS);