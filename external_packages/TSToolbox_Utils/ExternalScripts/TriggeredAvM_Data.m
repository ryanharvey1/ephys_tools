%function [EegSegAv EegSegStd EegSTrange] = TriggeredAvM(data,T,chanIx,win,sr,nChannels,SignalType)
% trigered averages subset of channels from file or matrix
% within a window of win (in msec)
% method : 1, to load all channels for each segment around trigger (if nChannels/nTriggers is large)
%               2, to load each channel and trigger average it separately (use that if nChannels/nTriggers is small)
%From A Sirota, changed by A Peyrache, 2012


function [EegSegAv, EegSegStd]=TriggeredAvM_Data(data,T,varargin)
[win,sr] = DefaultArgs(varargin,{1, 20000});

if (length(win)==1)
    win = [win win];
end
swin = round(win*sr/1000);
Srange = [-swin(1):swin(end)];
nsamples = length(Srange);
%Trange = linspace(-swin(1),swin(end),length(Srange))*1000/sr;

maxSample = length(data);
T = T(find(T>swin(1)+1 & T<maxSample-swin(end)));
EegSegAv = zeros(length(Srange),1);
EegSegStd = zeros(length(Srange),1);
warning off
[EegSegAv EegSegStd] = TriggeredAv(data, swin(1),swin(end),T);
warning on
