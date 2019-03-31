function samplerate=get_framerate(ts)
% get_framerate: calculates video frame rate
%
%   Input
%           ts: video timestamps in microseconds (neuralynx standard)
%   Output
%           samplerate: sample rate in Hz
%
% Ryan harvey 2019

tempts=(ts-ts(1))/10^6;
idx=find(tempts>=1);
samplerate=idx(1)-1;
end