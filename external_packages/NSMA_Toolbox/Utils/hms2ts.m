function ts1 = hms2ts(hh,mm,ss)

% hms2ts  Converts time in hours,minutes,seconds into timestamp (0.1 msec) units
%
% ts1 = hms2ts(hh,mm,ss)
%
% INPUTS:
%       hh - hours
%       mm - minutes
%       ss - seconds
% OUTPUTS:
%       ts1 - time in timestamp (0.1 msec) units


ts1 = (3600*hh + 60*mm + ss)*10000;
