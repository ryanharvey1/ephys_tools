function [ VelFiltered ] = VelFilter( ts,vel_cmPerSec,MinRatSpeed)
%VelFilter velocity filters data
%   Input: 
%           Matrix:         [TS,X,Y,etc. ]
%           vel_cmPerSec:   Instaneous velocity in cm/sec
%           MinRatSpeed:    Min velocity wanted
%   Output:
MinRunTime = 30;     % minimum time to be considered a running period (1 sec, in ts units)
MinPauseTime = 30;   % minimum time to be considered a break in running (1 sec, in ts units)

% X=matrix(:,2);
vel_t=ts(:,1);

% find times when rat is running slower than MinRatSpeed (cm/s)
istops=find(vel_cmPerSec<MinRatSpeed);
% stops_ts=vel_t(istops);

% find intervals when rat is not running
icross=find(diff(istops)>2);
start_i = [istops(1); istops(icross+1)];
end_i = [istops(icross); istops(end)];

S=[vel_t(start_i); vel_t(end)];
E=[vel_t(1); vel_t(end_i)];

% cut out runs shorter than MinRunTime
ii = 1;
while ii < length(S)
  while (S(ii) - E(ii)) < MinRunTime && ii < length(S)
    S(ii) = [];
    E(ii) = [];
  end
  ii = ii+1;
end

% cut out pause intervals shorter than MinPauseTime
ii=2;
while ii <= length(S)
  while (E(ii)-S(ii-1)) < MinPauseTime && ii < length(S)
    S(ii-1) = [];
    E(ii) = [];
  end
  ii = ii+1;
end

% figure; plot(vel_t, X, 'k-')
% xlabel('Time (ts from start of cheetah)')
% ylabel('Position (pixels)')
% axis tight
% hold on; plot(stops_ts, interp1(vel_t, X, stops_ts), 'r.')
% hold on; plot(S, interp1(vel_t, X, S), 'g.')
% hold on; plot(E, interp1(vel_t, X, E), 'b.')

startrun=E(1:end);
stoprun=S(1:end);

VelFiltered=[];
for i=1:length(startrun)
    VelFiltered=[VelFiltered;ts(find(ts(:,1)==startrun(i),1):find(ts(:,1)==stoprun(i),1),:)];
end
end

