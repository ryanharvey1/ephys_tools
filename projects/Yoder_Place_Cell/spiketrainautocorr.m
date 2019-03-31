% spiketrainautocorr
clear;clc;close all
addpath ('/Users/ryanharvey/GoogleDrive/MatlabDir/CircStat2012a','/Users/ryanharvey/GoogleDrive/MatlabDir/BClarkToolbox/Analysis')
load('/Users/ryanharvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/NewOCC_Map_workspace2.mat')
FigureLocation='/Users/RyanHarvey/Dropbox/school work/UNM/Lab/Projects/Place cells project with Yoder/Figures_For_Paper';
%%
close all
sessions=fieldnames(dataC);
sessions=sessions([1:5:length(sessions)],1);
i=32;
% for i=1:length(sessions)
frames=dataC.(sessions{i});
spk=frames(frames(:,5)==1,:);
% isi=diff(spk(:,1));
% R = corrcoef(isi,isi);
%%
ts=spk(:,1)/100000;

max_lag = 0.5;
t_bin=0.01;

% Acor - taken from intrinsic frequency 2
if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')
    max_lag = t_bin*floor(max_lag/t_bin)+.5*t_bin;
end

[cor, lag] = CrossCorr(ts, ts, 'lag', [-max_lag max_lag], 'binsize', t_bin, 'norm', 'count');


cor = cor/max(cor(lag~=0))*2-1;
S = abs(fft(cor)).^2;
df = 1/range(lag); fNQ = 1/mode(diff(lag))/2;
f = 0:df:fNQ;
S = S(1:length(f));
% The power spectrum was smoothed with a 2-Hz rectangular window
S = smooth(S,2/df);

% and the peak value in the 5-11 Hz band was identified
peak = f(f>=5&f<=11);
[~,i] = max(S(f>=5&f<=11));
peak = peak(i);

% A neuron was defined as theta)modulated if the mean power within 1)Hz of 
% each side of the peak in the 5�11 Hz frequency range was at least 5 times
% greater than the mean spectral power between 0 Hz and 50 Hz
ind = mean(S(abs(f-peak)<1))/mean(S(f<50));

%%
t1=spk(:,1)';
t2=spk(:,1)';
M1 = length(t1);
M2 = length(t2);
D = ones(M2,1)*t1 - t2'*ones(1,M1);
D = D(:);
D=D/1000000;
% figure;hist(D,200)

% end
[N,edges] = histcounts(D,-500:500);
N=smooth(N,20);
% figure;plot(-500:499,N)


%     opts = fitoptions('Method', 'NonlinearLeastSquares');
%     opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf 4];
%     opts.StartPoint = [0 0 0 0 0 0 0 7];
%     opts.Upper = [Inf Inf Inf Inf Inf Inf Inf 12];
%
% [f2, G] = fit([-500:499]',N,'fourier8',opts)
%  Frequency=f.w;
%     Fit=G.rsquare;

opts = fitoptions('fourier8');
%     opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf 4];
opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf 4];

%     opts.StartPoint = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7];
opts.Upper = [Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf 12];

[f2, G] = fit([-500:499]',N,'fourier8',opts)

Frequency=f2.w
Fit=G.rsquare
plot(f2,[-500:499]',N)
title(['Estimated frequency = ', num2str(Frequency), ' Hz']);
%%
opts = fitoptions('fourier1');
%     opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf 4];
opts.Lower = [-Inf -Inf -Inf 4];

opts.StartPoint = [0 0 0 7];
opts.Upper = [Inf Inf Inf 12];

[f2, G] = fit([-500:499]',N,'fourier1',opts)

Frequency=f2.w
Fit=G.rsquare
plot(f2,[-500:499]',N)
title(['Estimated frequency = ', num2str(Frequency), ' Hz']);


%%
t=[-500:499]';
x=N;
[F, G] = fit(t, x, 'fourier8')

% Plot result
plot(F,t,x);
grid on;
xlabel('Time (ms)');
ylabel('Amplitude');
title(['Estimated frequency = ', num2str(F.a8), ' Hz']);
%%


Fs = 1000;                   % samples per second
dt = 1/Fs;                   % seconds per sample
StopTime = 1;             % seconds
t = (0:dt:StopTime-dt)';     % seconds
% Sine wave:
Fc = 1:100;                     % hertz
%    x = cos(2*pi*Fc*t);
x=cos(2*pi*t*Fc);

% glmCoef=zeros(9);
% numModel = 9;
% for iModel = 1 : numModel
% clear curPredSpace;
% curPredSpace= x'; %create a temporary predictor matrix
% curPredSpace([iModel:numModel],:) = []; % remove the predictors for the current model from the predictor matrix;
% [glmCoef(iModel,:)] = glmfit(curPredSpace', N)
[b,dev,stats]=glmfit(x, N)

% end



% t=N;
% fun=@(t)[a*(sin(wt)+1)+b]*(e^-abs(t)/T1)+c*e^-t^2/T2^2;
%
% x = fminsearch(fun)

%%
% [autocor,lags] = xcorr(diff(spk(:,1)));
% figure;plot(lags/60,autocor)
% xlabel('Lag (days)')
% ylabel('Autocorrelation')
binsize = 4;
nbins = 250;
[h, t] = AutoCorr((spk(:,1)/1000000), binsize, nbins);
hi = h./max(h);
figure (1), plot(t, hi, 'k'),
ylabel('Normalized Auto-Correlation')
xlabel('Lag (msec)')
xlim([0 500]);
box off
%%
%load spike train data from .ts file outputs (.ts.r)
ts = spk(:,1); %load spike data

%pad ts to have enough cells to run iterative code
pad=zeros(150,1);
ts_pad=[ts;pad];

%convert spikes from 10s of usec to msec
ts_msec = ts_pad./1000000;

ISI_list=[];
ISIx=[];
ISIy=[];
y=([2:101])';
%subtract each ts from all subsequent ts to generate ISI list
for j = 1:(length(ts_pad)-149);
    ISIx = ts_msec(j+1)-ts_msec(j);
    ISI_list = [ISI_list; ISIx];
    if ISIx >= 0 && ISIx <= 501;
        ISIy = ts_msec(j+y)-ts_msec(j);
        ISI_list = [ISI_list; ISIy];
    else ISIx >= 500;
    end
end

%remove all negative values, zeros, and values over 500
ISIfinal0=ISI_list;
ISIfinal = ISIfinal0(ISIfinal0~=0);
ISIfinal(ISIfinal <= 0) = [];
ISIfinal(ISIfinal > 500) = [];

%bin ISI by delay
[ISIbins]=hist(ISIfinal,(1:500));

%smooth autocorrelation curve and normalize to 1
ISIsmooth = smooth(ISIbins);
ISIsmooth(ISIsmooth < 0) = [0];
maxISIsmooth = max(ISIsmooth);
normISIsmooth = ISIsmooth./maxISIsmooth;

%calculate a theta ratio (first trough/second peak)
minISI = min(normISIsmooth(25:105));
maxISI = max(normISIsmooth(95:175));
ThetaRatio = (maxISI-minISI);

%fit a sinusoid to the autocorrelation
xSinu=[0.00628318:0.00628318:3.14159]';
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Lower = [-Inf -Inf -Inf 4];
opts.StartPoint = [0 0 0 7];
opts.Upper = [Inf Inf Inf 12];
[f, G]=fit(xSinu,normISIsmooth,'fourier1',opts);
Frequency=f.w;
Fit=G.rsquare;

%plot the autocorrelation
fig1 = area(normISIsmooth,'facecolor','k'); %smoothed and normalized autocorrelation
axis square tight
set(gca,'XTick',[0 100 200 300 400 500],'YTick',[0 0.5 1.0]);
box off
title(['Frequency = ', num2str(Frequency), ' Hz', ' | Fit= ',num2str(Fit)]);
hold on
plot(f)

function [cor, lags] = CrossCorr(ts1, varargin)
% [acor, lag] = CMBHOME.Spike.CrossCorr(ts1)
% [xcor, lag] = CMBHOME.Spike.CrossCorr(ts1, ts2, varargin)
%
% Calculates cross correlation (or auto correlation) for vectors ts1 and
% ts2 of spike times. ts1 and ts2 are not binned, or binary, spike trains but
% rather the time that spikes occured.
%
% If ts1 is the only argument, the autocorrelation is produced (also if
% ts1==ts2). In this case, all zeros in the set of latencies will be removed
% before calculating the xcorr histogram.
%
% ARGUMENTS
%   ts1         vector of spike times (seconds)
%   ts2         (optional) vector of spike times
%
% PROPERTIES
%   'binsize'       (default=.01 seconds) binsize for cross correlation
%   'norm'          (default='prob') determines normalization of cross
%                   correlogram. Can be 'prob', 'count', 'unbiased'. To use
%                   the unbiased option, 'epoch' property must be passed.
%   'lag'           (default [-.4 4] seconds) included to speed up algorithm.
%                   this defines the upper and lower limit in the peristumulus
%                   time histogram
%   'suppress_plot' (default 1) 1 with plot, 0 without plot
%   'epoch'         2 element vector indicating start and stop times of recording window
%
% RETURNS
%   cor         a col. vector of cross correlation values corresponding to 'lag'
%   lag         a col. vector of center-aligned lag values (seconds)
%
% alex wiltschko and andrew bogaard
% Updated 3rd August 2012, Ehren Newman, increased speed for large
% sessions.

p = inputParser;

p.addRequired('ts1');
p.addOptional('ts2', [], @(c) isnumeric(c));
p.addParamValue('norm', 'prob', @(c) ischar(c));
p.addParamValue('binsize', .01, @(c) numel(c)==1 && (c>0));
p.addParamValue('lag', [-.4 .4], @(c) numel(c)==2 && diff(c)>0);
p.addParamValue('suppress_plot', 1, @(c) numel(c)==1);
p.addParamValue('epoch', [], @(c) numel(c)==2);

p.parse(ts1, varargin{:});

ts1 = p.Results.ts1(:); % make col vectors
ts2 = p.Results.ts2(:);
norm = p.Results.norm;
binsize = p.Results.binsize;
lag = p.Results.lag;
suppress_plot = p.Results.suppress_plot;
epoch = p.Results.epoch;

ac = 0;

if isempty(ts2), ts2 = ts1; end

if isequal(ts1,ts2), ac = 1; end % is autocorr

db = nan(length(ts1), 3);

s1 = 1;

spkind = 1;

while spkind <= length(ts1)
   
    s = s1;
    
    while ts2(s) < ts1(spkind)+lag(1) && s < length(ts2)
    
        s = s+1;
    
    end
    
    s1 = s;
    
%    f = s;
        
    f = find(ts2 <= ts1(spkind)+lag(2)==0,1);
    
    if isempty(f), f = length(ts2)+1; end
      
%     while ts2(f) <= ts1(spkind)+lag(2)
%         
%         f = f+1;
%         
%         if f>length(ts2), break; end
%         
%     end
        
    if ts2(s)<=ts1(spkind)+lag(2)
        db(spkind, :) = [s f-1 ts1(spkind)];
    end
    
    spkind = spkind+1;
    
end

dspk = diff(db(:,1:2), 1, 2);

N = 0;
for i = 0:max(dspk)
    
    N = N + sum(dspk>=i);
    
end


lags = lag(1)+binsize/2:binsize:lag(2)-binsize/2;

saveMem = N*4 > 1000000000;
if saveMem
  lags = sort([lags 0-1e-10 0 0+1e-10]);
  cor = single(nan(max(dspk)+1,length(lags)));
else
  psth = single(ones(N,1));
end

N = 1;
for i = 0:max(dspk)
    
    where = dspk>=i;
    
    tf = db(where,1)+i;
    
    if saveMem cor(i+1,:) = single(hist(ts2(tf)-db(where,3),lags));
    else psth(N:N+sum(where)-1) = single(ts2(tf)-db(where,3)); end

    N = N+sum(where);
    
end

if saveMem
  cor = sum(cor);
  zeroLag = find(lags==0);
  cor = [cor(1:zeroLag-3),sum(cor(zeroLag-2:zeroLag-1),2), sum(cor(zeroLag+2:zeroLag+1),2), cor(zeroLag+3:end)];
else
  if ac, psth(psth==0) = []; end % remove zeros in autocorrelation
  psth(psth<lags(1)-binsize/2) = [];
  cor = hist(psth, lags);
end

if strcmp(norm, 'unbiased')
   
    if isempty(epoch)
        
        disp('''epoch'' must be defined for unbiased normalization, no normalization performed')
    
    else
        
        L = floor(diff(epoch)/binsize);

        normc = abs(lags)/binsize;
        
        cor = cor./(L-normc);
        
    end
    
elseif strcmp(norm, 'prob')
    
    cor = cor / min(length(ts1), length(ts2));
    
end % otherwise, it will be count

if ~suppress_plot
    bar(lags*1000, cor, 1, 'k'), xlim(1000*[lag(1) lag(2)]), hold on
    line([0 0], ylim, 'linestyle', ':', 'color', [.7 .7 .7]);
    set(gca, 'fontsize', 8, 'box', 'off');
end

cor = cor(:);

lags = lags(:);
end

