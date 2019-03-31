function [hh, mm, ss, hms_string] = ts2hms(ts1)

% ts2hms  Converts time in timestamp (0.1 msec) units into hours,minutes,seconds
%
% [hh, mm, ss, hms_string] = ts2hms(ts1)
%
% INPUTS:
%       ts1 - time in timestamp (0.1 msec) units
% OUTPUTS:
%       hh - hours
%       mm - minutes
%       ss - seconds
%       hms_string - time in string format hh:mm:ss.ssss


hh = floor(ts1/(3600*10000));
mm = floor(mod(ts1,3600*10000)/(60*10000));
ss = mod(ts1,60*10000)/10000;
hms_string = [num2str(hh) ':' num2str(mm) ':' num2str(ss)];
