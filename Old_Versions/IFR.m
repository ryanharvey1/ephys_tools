function ifr=IFR(data_video,spikets,samplerate)
% IFR -instantaneous firing rate 
% 
% INPUT:
%           spikets:    timestamps for each spike
%           spk_binary: binary representing each spike
%           sessdur:    duration of session on seconds
%           samplerate: sampling rate 
% 
% OUTPUT:
%           ifr: instantaneous firing rate
%
% RYAN HARVEY 2017, updated 2019

% seconds to bin over
binsize=.2;
% bin over .200 seconds
edges=0:binsize:data_video(end,1);
[N,~] = histcounts(spikets,edges);
% convert to spikes/sec
N=N*(samplerate/round(binsize*samplerate));
% smooth 
ifr=smooth(N,round(binsize*1.5*samplerate));
end

