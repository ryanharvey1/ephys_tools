function Sres = RestrictTheSMatrix2Interval(S,offset_sec,interval_sec)
% Restrict a S-Matrix (Spike Matrix) S to an interval starting at offset_sec with duration of length interval_sec.
%
% Sres = RestrictTheSMatrix2Interval(S,offset_sec,interval_sec)
% 
% This function is intended to be used in spleep analysis where one typically
% restricts all cells in sleep1 to the last 10 min of sleep and in sleep2 to the 
% first 10 mins after the behaviour epoch.
% The offset_sec allows to shift the interval a few minutes forward (positve offset) or backwards
% (negative offset, counted from the end).
% Both offset_sec and interval_sec are given in seconds.
%
% First, the beginning  and end  times of the session (SessStartTS, SessStopTS) are determined by 
% finding the mins and maxs of all timestamps in all cells (ts-objects) in S.
%
% If interval_sec is positive, the interval is taken from the start of the session after abs(offset_sec):
% => Sres is then restricted to [SessStartTS+abs(offset_sec)*10000, SessStartTS+(abs(offset_sec)+abs(interval_sec))*10000];
%
% If interval_sec is negative, the interval is taken from the end of the session and shifted forward by offset_sec:
% => Sres is then restricted to [SessEndTS-(abs(offset_sec)+abs(interval_sec))*10000, SessEndTS-abs(offset_sec)*10000];
%
% PL Feb 2003

% find SessStartTS, SessStopTS:
[SessStartTS, SessEndTS] = SessionStartEndTS(S);

% determine restriction interval [ts0,ts1]: 
if interval_sec > 0
    ts0 = SessStartTS + abs(offset_sec)*10000;
    ts1 = SessStartTS + (abs(offset_sec)+abs(interval_sec))*10000;
else
    ts0 = SessEndTS - (abs(offset_sec)+abs(interval_sec))*10000;
    ts1 = SessEndTS - abs(offset_sec)*10000;
end

Sres = RestrictTheSMatrix(S, ts0, ts1);