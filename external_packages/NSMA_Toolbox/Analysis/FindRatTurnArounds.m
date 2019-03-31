function [startgoodrun, stopgoodrun]=FindRatTurnArounds(V,laps,runthresh)
% find times when rat turns around in the middle of a lap
% inputs: V: coordinates of rat, parametrized to track
%         laps: struct with begin times and directions of each lap
%         runthresh: threshold velocity considered running 
%                    (units same as V / sec)
%                    any velocities slower than this will not be considered
% outputs: startgoodrun, stopgoodrun: start and stop times of periods of 
%          running in correct direction
% Running in good direction periods are greater than 0.5s, and running in 
% bad direction periods are longer than 1s
% 
% ZN 04/2011

MinGoodDirectionTime = 5000;     % minimum time to be considered a good running period (1 sec, in ts units)
MinBadDirectionTime = 10000;   % minimum time to be considered a break in running (1 sec, in ts units)
framesPerSec=30;            % tracker sampling rate

% smooth position data
sV=SmoothTsd(V,framesPerSec);
% find velocity of rat from track position
vel_v = diff(Data(sV));			% vel units are pixels/frame
vel_t = Range(V,'ts');
vel_t(1) = [];				% need to delete first point to make vel_t same length as vel_v 
%v=Data(V);
%v(1)=[];

% find times when rat is running in wrong direction of the lap
%a=figure;
wrongdir_idx=[];
vel_idx=[1:length(vel_t)]';
for l=1:length(laps)
    startlap=laps(l).start_ts;
    if l<length(laps)
        endlap=laps(l+1).start_ts;
    end

    vlap = vel_v(vel_t<endlap & vel_t>startlap);
    vlapdir = sign(vlap);
    idxlap = vel_idx(vel_t<endlap & vel_t>startlap);
    %figure(a); plot(idxlap, v(vel_t<endlap & vel_t>startlap), 'k-')   
    baddiridxlap = vlapdir~=laps(l).direction & abs(vlap*framesPerSec)>runthresh;
    %hold on; plot(idxlap(baddiridxlap), v(idxlap(baddiridxlap)), 'r.')
    wrongdir_idx = [wrongdir_idx; idxlap(baddiridxlap)];
    %title(['lap # ', num2str(l), ' direction: ', num2str(laps(l).direction)])
    %pause;
    %hold off;
end
wrongdir_t = vel_t(wrongdir_idx);

% find intervals when rat is turned around
iturn = find(diff(wrongdir_idx)>2);
start_i = [wrongdir_idx(1); wrongdir_idx(iturn+1)];
end_i = [wrongdir_idx(iturn); wrongdir_idx(end)];

S=vel_t(start_i);
E=vel_t(end_i);

% cut out gaps shorter than MinGoodDirectionTime
ii = 2;
while ii < length(S)
  while (S(ii) - E(ii-1)) < MinGoodDirectionTime && ii < length(S)
    S(ii) = [];
    E(ii-1) = [];
  end
  ii = ii+1;
end

% cut out theta intervals shorter than MinBadDirectionTime
ii=1;
while ii <= length(S)
  while (E(ii)-S(ii)) < MinBadDirectionTime && ii < length(S)
    S(ii) = [];
    E(ii) = [];
  end
  ii = ii+1;
end

figure; plot(Range(V), Data(V), 'k-')
xlabel('Time (ts from start of cheetah)')
ylabel('Position on Track')
axis tight
hold on; plot(wrongdir_t, interp1(Range(V), Data(V), wrongdir_t), 'r.')
hold on; plot(S, interp1(Range(V), Data(V), S), 'g.')
hold on; plot(E, interp1(Range(V), Data(V), E), 'b.')

startgoodrun=[vel_t(1); E(1:end)];
stopgoodrun=[S(1:end); vel_t(end)];
