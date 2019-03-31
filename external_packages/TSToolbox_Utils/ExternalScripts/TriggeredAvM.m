%function [EegSegAv EegSegStd EegSTrange] = TriggeredAvM(FileBase,T,win,sr,nChannels,method, SignalType)
% trigered averages multiple channels from file or matrix
% within a window of win (in msec)
% method : 1, to load all channels for each segment around trigger (if nChannels/nTriggers is large)
%               2, to load each channel and trigger average it separately (use that if nChannels/nTriggers is small)

function [EegSegAv, EegSegStd, Trange]=TriggeredAvM(Filebase,T,varargin)
[win,sr,nChannels,method, SignalType] = DefaultArgs(varargin,{1000, 1250, [], 1, 'eeg'});

if (length(win)==1)
    win = [win win];
end
swin = round(win*sr/1000);
Srange = [-swin(1):swin(end)];
nsamples = length(Srange);
Trange = linspace(-swin(1),swin(end),length(Srange))*1000/sr;


if isempty(nChannels)
    Par = LoadPar(Filebase);
    nChannels =Par.nChannels;
end

maxSample = FileLength([Filebase '.' SignalType]) /nChannels/2;
T = T(find(T>swin(1)+1 & T<maxSample-swin(end)));
EegSegAv = zeros(length(Srange),nChannels);
EegSegStd = zeros(length(Srange),nChannels);
if method == 1
    nT = length(T);
    for d=1:nT
        startpos = (T(d)-swin(1))*nChannels*2;
        tmpeeg = bload([Filebase '.' SignalType], [nChannels nsamples], startpos,'int16');
        tmpeeglong = bload([Filebase '.' SignalType], [nChannels nsamples*5], startpos,'int16');
        tmpeeg = tmpeeg - repmat(mean(tmpeeglong,2),1,nsamples);
        EegSegAv = EegSegAv+tmpeeg(:,:)';
        EegSegStd = EegSegStd+tmpeeg(:,:)'.^2;
    end
    EegSegAv = EegSegAv./nT;
    EegSegStd = EegSegStd./nT - EegSegAv.^2;
else

    for i=1:nChannels
        eeg = LoadBinary([Filebase '.' SignalType], i, nChannels,2)';
        [EegSegAv(:,i) EegSegStd(:,i)] = TriggeredAv(eeg, swin(1),swin(end),T);
    end
end
