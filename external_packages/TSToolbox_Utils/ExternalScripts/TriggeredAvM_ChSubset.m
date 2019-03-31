%function [EegSegAv EegSegStd EegSTrange] = TriggeredAvM(data,T,chanIx,win,sr,nChannels,SignalType)
% trigered averages subset of channels from file or matrix
% within a window of win (in msec)
% method : 1, to load all channels for each segment around trigger (if nChannels/nTriggers is large)
%               2, to load each channel and trigger average it separately (use that if nChannels/nTriggers is small)
%From A Sirota, changed by A Peyrache, 2012


function [EegSegAv, EegSegStd, Trange]=TriggeredAvM_ChSubset(Filebase,T,chanIx,varargin)
[win,sr,nChannels,SignalType] = DefaultArgs(varargin,{1000, 20000, [], 'dat'});

if (length(win)==1)
    win = [win win];
end
swin = round(win*sr/1000);
Srange = [-swin(1):swin(end)];
nsamples = length(Srange);
Trange = linspace(-swin(1),swin(end),length(Srange))*1000/sr;

if isempty(nChannels)
    Par = LoadXml(Filebase);
    nChannels =Par.nChannels;
end

maxSample = FileLength([Filebase '.' SignalType]) /nChannels/2;
T = T(find(T>swin(1)+1 & T<maxSample-swin(end)));
EegSegAv = zeros(length(Srange),length( ));
EegSegStd = zeros(length(Srange),nChannels);

chanIx = chanIx(:)';

for ii=chanIx
    eeg = LoadBinary([Filebase '.' SignalType],'channels',ii,'nChannels',nChannels,'frequency',sr)';
    [EegSegAv(:,ii) EegSegStd(:,ii)] = TriggeredAv(eeg, swin(1),swin(end),T);
end
