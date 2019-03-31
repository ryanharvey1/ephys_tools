function [startrun, stoprun, Vel_cmps, Vel_phi]=FindRatRunsXY(X,Y,pixels2cm,MinRatSpeed,framesPerSec)
% [startrun, stoprun, Vel_cmps, Vel_phi]=FindRatRunsXY(X,Y,pixels2cm,MinRatSpeed,framesPerSec)
% 
% find times when rat is running in any environment (not stopped)
% inputs: X,Y: coordinates of rat
%         pixels2cm: conversion from pixels to cm
%         MinRatSpeed: minimum speed of rat (in cm/s) to consider running
%         framesPerSec: video tracker sampling rate ~30 fps for digital
%                       cheetah recordings, ~60 fps for analog cheetah
%                       default = 30 fps
% outputs: startrun, stoprun: start and stop times of running periods
%          Vel_cmps: velocity of rat (cm/s)
%          Vel_phi: angle phi of velocity vector = "direction of velocity" in radians   
% Running periods are greater than 2s, and stopped periods are longer than 1s
% 
% ZN 2010 

% (used parts from Peter's code (from his course, and FindThetaStates))

MinRunTime = 20000;     % minimum time to be considered a running period (1 sec, in ts units)
MinPauseTime = 10000;   % minimum time to be considered a break in running (1 sec, in ts units)
if nargin<5
    framesPerSec = 30;      % tracker sampling rate  
end

% smooth position tsd's
Xs = SmoothTsd(X, framesPerSec);		% smooth X with 1 sec smoothing window
Ys = SmoothTsd(Y, framesPerSec);		% smooth Y with 1 sec smoothing window
Xs=Restrict(Xs, StartTime(X)+10000, EndTime(X)-10000);     % don't use first and last second, smoothing doesn't work here
Ys=Restrict(Ys, StartTime(X)+10000, EndTime(X)-10000);

figure; plot(Data(Xs), Data(Ys), 'k-')
title('Smoothed Positions')
xlabel('Xs')
ylabel('Ys')
axis equal

% find velocity of rat from smoothed position
vel_x = diff(Data(Xs));			% vel units are pixels/frame
vel_y = diff(Data(Ys));
vel_t = Range(Xs,'ts');
vel_t(1) = [];				% need to delete first point to make vel_t same length as vel_x and vel_y  

% figure; 
% plot(vel_t/10000/60, [vel_x, vel_y])
vel_abs = sqrt( vel_x.^2 + vel_y.^2);		%  scalar length of velocity vector = "scalar velocity" in pixels/frame
vel_phi = atan2(vel_y, vel_x);			%  angle phi of velocity vector = "direction of velocity" in radians   

vel_cmPerSec = vel_abs * pixels2cm * framesPerSec;  

% figure;
% plot(vel_t/60000, vel_cmPerSec)

% find times when rat is running slower than MinRatSpeed (cm/s)
istops=find(vel_cmPerSec<MinRatSpeed);
stops_ts=vel_t(istops);

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

figure; plot(Range(X), Data(X), 'k-')
xlabel('Time (ts from start of cheetah)')
ylabel('Position (pixels)')
axis tight
hold on; plot(stops_ts, interp1(Range(X), Data(X), stops_ts), 'r.')
hold on; plot(S, interp1(Range(X), Data(X), S), 'g.')
hold on; plot(E, interp1(Range(X), Data(X), E), 'b.')

startrun=E(1:end);
stoprun=S(1:end);

Vel_cmps=tsd(vel_t, vel_cmPerSec);
Vel_phi=tsd(vel_t, vel_phi);