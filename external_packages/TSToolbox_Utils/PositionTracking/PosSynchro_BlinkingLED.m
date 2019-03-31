function PosSynchro_BlinkingLED(fbasename,nchannels,LedChannel,varargin)

if ~isempty(varargin)
    colorvec = varargin{1};
else
    colorvec = [1 3];
end

syncS = LoadBinary([fbasename '.dat'],'frequency',20000,'nchannels',nchannels,'channels',LedChannel);

tDat = (0:length(syncS)-1)/2;
epDat = thresholdIntervals(tsd(tDat,syncS),1000);
stDat = Start(epDat);
tFrstDat = stDat(1);
tLastDat = stDat(end);

tsp = load([fbasename '.tsp']);
tTsp = (0:size(tsp,1)-1);
epTsp = thresholdIntervals(tsd(tTsp,tsp(:,4)),500);
epTsp = mergeCloseIntervals(epTsp,5);
epTsp = dropShortIntervals(epTsp,5);

stTsp = Start(epTsp);

if length(stTsp) ~= length(stDat);
    keyboard
end

tFrstTsp = stTsp(1);
tLastTsp = stTsp(end);

totDur = tLastDat - tFrstDat;
videoSamplingRate = (tLastTsp-tFrstTsp) / (tLastDat-tFrstDat);

tTsp    = tTsp / videoSamplingRate;
stTsp   = stTsp / videoSamplingRate;
dt = median(stTsp-stDat);

tTsp = tTsp - dt;
stTsp = stTsp- dt;

if max(stDat-stTsp)>350
    keyboard
end

t = (0:256:tDat(end));
interpTsp = zeros(length(t),7);

interpTsp(:,1)=t;
warning off
    tsp(tsp==-1)        =   NaN;
    interpTsp(:,2:end)  =   interp1(tTsp,tsp(:,2:end),interpTsp(:,1));
warning on
interpTsp(isnan(interpTsp))=-1;


if 1
    figure(1),clf
    scatter(interpTsp(:,2),interpTsp(:,3),'r.')
    hold on
    scatter(interpTsp(:,6),interpTsp(:,7),'b.')
    %pause
    %keyboard
end

colorvec = sort([2*colorvec-1 2*colorvec])+1;
dlmwrite([fbasename '.whl'],interpTsp(:,colorvec),'delimiter','\t');