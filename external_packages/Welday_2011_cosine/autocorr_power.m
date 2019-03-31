% This script can be used to analyze the example datafile:
%
% Fig1_rawdata.mat - contains raw spike and tracking data for the example
%                    cell in Fig. 1 of Welday et al.
%
% It may also be used to analyze properly formatted data supplied by the
% user.
% -------------------------------------------------------------------------
%
% Before executing this script, the following variables must be assigned,
% either by loading in the example data file or by loading in user 
% data that adheres to the following format:
%
% Position_Speed(N,2) -- each of the N rows is a position sample, and
% the two matrix columns are:
%   Position_Speed(:,1)=running speed in units of pixels per sample
%   Position_Speed(:,2)=time stamp in units of seconds
%
% spikedata -- theta cell spike time stamps (in seconds)
% lobound -- timestamp (in seconds) at which data analysis begins
% hibound -- timestamp (in seconds) at which data analysis ends
%
% e_fints, ne_fints, n_fints, nw_fints, w_fints, sw_fints, s_fints, se_fints
% These are two-column lists of movement intervals (first column start,
% second column end) for each direction, timestamped in the same base as
% the speed data in 'Position_Speed'
%
% -------------------------------------------------------------------------
% 
% In addition to the above variables, the example data file also contains 
% a variable called Position_Xcol1_Ycol2, in which the first and second
% columns are the X and Y pixel coordinates, respectively, at each time step. 
% The timestamps for each position sample may be found in the second column
% of Position_Speed. Although the X and Y pixel coordinates are not
% required for the analysis below, they are included for completeness of
% the example dataset.

% -------------------------------------------------------------------------
% The program generates output in six figure windows:
% Figure 1:  The DBFT curve of the cell 
% Figure 19: The spike autocorrelogram for all movement directions averaged together
% Figure 20: The burst frequency power spectrum for each direction, as in Fig. 1E of Welday et al.
% Figure 21: The spike autocorrelograms for each direction, as in Fig. 1D of Welday et al.
% Figure 32: The running speed distributions for each direction, with balanced distribution in the
%            center, as in Fig. 1B of Welday et al.
% Figure 34: The spike rate histograms (25 m bins) for movement in each direction (all movement
%            epochs for each direction are concatenated together in sequential order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%ASSIGN ANALYSIS PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

samplerate = 30; %position sample rateFi
pixpercm = 4.7; %pixel per cm
blocksec=.4; %%size of movement blocks to analyze (in seconds), during which direction must be constant & speed must meet criterion
blocksamp=round(blocksec/(1/samplerate)); %%position samples per blocksec
binsize=.025; %firing rate histogram bin width (in seconds)

%accumulate results in these variables:
SpeedDists=[]; %each row is a distribution of running speeds (indexed by 'speed_edges') for a different direction
speedbydir=[]; %each element is the expected running speed in a different direction, to summarixe how running speed varies with direction
BalSpeedDists=[];
meanspeeds=[];
speedstats=[];

%%%restrict analysis to specified time boundaries
Position_Speed=Position_Speed(find((Position_Speed(:,2)>=lobound) & (Position_Speed(:,2)<=hibound)),:);

%%edges and centers of bins for computing running speed distributions
%%(specified in units of cm/s)
speed_edges=[5 7.5 10 12.5 15 17.5 20 22.5 25 27.5 30 32.5 35 37.5 40 42.5 45 47.5 50];
speed_middles=[5 7.5 10 12.5 15 17.5 20 22.5 25 27.5 30 32.5 35 37.5 40 42.5 45 47.5]+1.25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s%%%%%%%%%%%%%%%
%%%%%COMPUTE SPEED DISTRIBUTIONS & RATE HISTOGRAMS IN EACH DIRECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x, e_speeds, e_ratehisto, e_autohisto]=block_speeds(e_fints, Position_Speed, spikedata, pixpercm, blocksamp, binsize, speed_edges, blocksec, samplerate);
SpeedDists=[SpeedDists x];  %%accumulate the distribution of running speeds for this direction as a new row in SpeedDists

[x, ne_speeds, ne_ratehisto, ne_autohisto]=block_speeds(ne_fints,Position_Speed,spikedata,pixpercm,blocksamp,binsize,speed_edges, blocksec, samplerate);
SpeedDists=[SpeedDists x];  %%accumulate the distribution of running speeds for this direction as a new row in SpeedDists

[x, n_speeds, n_ratehisto, n_autohisto]=block_speeds(n_fints,Position_Speed,spikedata,pixpercm,blocksamp,binsize,speed_edges, blocksec, samplerate);
SpeedDists=[SpeedDists x];  %%accumulate the distribution of running speeds for this direction as a new row in SpeedDists

[x, nw_speeds, nw_ratehisto, nw_autohisto]=block_speeds(nw_fints,Position_Speed,spikedata,pixpercm,blocksamp,binsize,speed_edges, blocksec, samplerate);
SpeedDists=[SpeedDists x];  %%accumulate the distribution of running speeds for this direction as a new row in SpeedDists

[x, w_speeds, w_ratehisto, w_autohisto]=block_speeds(w_fints,Position_Speed,spikedata,pixpercm,blocksamp,binsize,speed_edges, blocksec, samplerate);
SpeedDists=[SpeedDists x];  %%accumulate the distribution of running speeds for this direction as a new row in SpeedDists

[x, sw_speeds, sw_ratehisto, sw_autohisto]=block_speeds(sw_fints,Position_Speed,spikedata,pixpercm,blocksamp,binsize,speed_edges, blocksec, samplerate);
SpeedDists=[SpeedDists x];  %%accumulate the distribution of running speeds for this direction as a new row in SpeedDists

[x, s_speeds, s_ratehisto, s_autohisto]=block_speeds(s_fints,Position_Speed,spikedata,pixpercm,blocksamp,binsize,speed_edges, blocksec, samplerate);
SpeedDists=[SpeedDists x];  %%accumulate the distribution of running speeds for this direction as a new row in SpeedDists

[x, se_speeds, se_ratehisto, se_autohisto]=block_speeds(se_fints,Position_Speed,spikedata,pixpercm,blocksamp,binsize,speed_edges, blocksec, samplerate);
SpeedDists=[SpeedDists x];  %%accumulate the distribution of running speeds for this direction as a new row in SpeedDists

%%%-----------------------------------
%%%      plot speed distributions
%%%-----------------------------------
ymax=max(max(SpeedDists)); %scale for y-axis
figure(32); clf; %plot the speed distributions in this figure
for i=1:8 %loop through the movement directions
    subplot(3,3,i+1*(i>4)); %arrange graphs geometrically according to direction, leave middle blank for balanced distribution
    bar(speed_middles,SpeedDists(1:length(speed_middles),i)*blocksec,1,'r'); %plot distribution bar graph
    set (gca,'XLim',[0 50],'YLim',[0 ymax]); axis square;  %set axis scales
end
balcounts=(min(SpeedDists')>0).*max(SpeedDists'); %compute the balanced speed distribution
subplot(3,3,5); bar(speed_middles,balcounts(1:length(speed_middles))*blocksec,1); %plot balanced speed distribution in the center
set (gca,'XLim',[0 50],'YLim',[0 ymax]);  axis square; %set axi scales

balexpect=0; %expectation value for balanced speed distribution (needed to compute predicted grid size)
for i=find(balcounts>0)
    balexpect=balexpect+speed_middles(i)*balcounts(i)/sum(balcounts);
end

%%%-----------------------------------
%%%      plot rate histograms
%%%-----------------------------------
figure(34); clf;  %plot the rate histograms in this figure
hh=e_ratehisto'; subplot(8,1,1); bar(hh(:),1); axis tight;
hh=ne_ratehisto'; subplot(8,1,2); bar(hh(:),1); axis tight;
hh=n_ratehisto'; subplot(8,1,3); bar(hh(:),1); axis tight;
hh=nw_ratehisto'; subplot(8,1,4); bar(hh(:),1); axis tight;
hh=w_ratehisto'; subplot(8,1,5); bar(hh(:),1); axis tight;
hh=sw_ratehisto'; subplot(8,1,6); bar(hh(:),1); axis tight;
hh=s_ratehisto'; subplot(8,1,7); bar(hh(:),1); axis tight;
hh=se_ratehisto'; subplot(8,1,8); bar(hh(:),1); axis tight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%PERFORM SPECTRAL ANALYSIS OF THETA BURST FREQ IN EACH DIRECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%-----------------------------
%% spectral analysis parameters
%%-----------------------------

figbase=0;  %base for figure numbering
x=0:(2*pi/8):2*pi;  %analysis uses 8 directional bins
maxfrequency = 20;  %power spectra will be computed from 0 - maxfrequency
numautobins = 256;  %number of bins on either side of the zero line in the autocorrelogram
peakwidth = 1.5; %bandwidth (in Hz) on either side of the theta peak across which to intgrate for computing expected frequency value
f = (numautobins/blocksec)*(0:(2^18))/2^19;  %frequency bins of the power spectrum

%%-----------------------------
%%  compute power spectra
%%-----------------------------
global iterat

iterat=[];

%%accumulate the power spectra in each direction
nex2=[]; %variable for accumulating power spectra in each direction
nex2=[nex2; autopowspect(maxfrequency,numautobins,blocksec,e_autohisto,e_speeds(:,5),balcounts)];
nex2=[nex2; autopowspect(maxfrequency,numautobins,blocksec,ne_autohisto,ne_speeds(:,5),balcounts)];
nex2=[nex2; autopowspect(maxfrequency,numautobins,blocksec,n_autohisto,n_speeds(:,5),balcounts)];
nex2=[nex2; autopowspect(maxfrequency,numautobins,blocksec,nw_autohisto,nw_speeds(:,5),balcounts)];
nex2=[nex2; autopowspect(maxfrequency,numautobins,blocksec,w_autohisto,w_speeds(:,5),balcounts)];
nex2=[nex2; autopowspect(maxfrequency,numautobins,blocksec,sw_autohisto,sw_speeds(:,5),balcounts)];
nex2=[nex2; autopowspect(maxfrequency,numautobins,blocksec,s_autohisto,s_speeds(:,5),balcounts)];
nex2=[nex2; autopowspect(maxfrequency,numautobins,blocksec,se_autohisto,se_speeds(:,5),balcounts)];
nex2=nex2'; %transpose


%%normalize and smooth the power spectra
for i=1:8 %loop through the directions
    nex2(:,i)=nex2(:,i)-min(nex2(:,i)); %shift the power spectrum to a zero baseline
    nex2(:,i)=smooth(nex2(:,i),14); %smooth the power spectrum with a 14-bin boxcar window
end

%%--------------------------------------------------
%%  extract burst frequencies from the power spectra
%%--------------------------------------------------

rangedex=find(f<11); %we will seek the peak of the power spectrum in a band below 11 Hz
rangedex=find(f(rangedex)>5); %and above 5 Hz
[maxpow, maxdex]=max(nex2(rangedex,:)); %find the amplitude and location of the theta peak
maxdex=maxdex+rangedex(1)-1; %convert the index of the peak's location into the same integer base as the 'f' axis
nexchopped=nex2; %make a copy of the power spectra to plot thresholding
nextops=nex2*0;  %make a blank power spectrum (all zeros) to plot thresholding

for i=1:8 %loop through the directions
    pwid=round(peakwidth/(f(2)-f(1))); %number of frequency bins over which to integrate for computing expected frequency value
    thresh=(maxpow(i)*.5); %set power spectrum threshold at half the peak amplitude
    peakdex=find(nex2((maxdex(i)-pwid):(maxdex(i)+pwid),i)>thresh)+maxdex(i)-pwid-1;  %indices of frequency bins over which to integrate for computing expected frequency value
    nexchopped(peakdex,i)=thresh; %make a copy of the power spectrum with values above the threshold 'chopped off' (for area graph plotting)
    nextops(peakdex,i)=nex2(peakdex,i)-thresh; %make another copy with only the above-threshold part (also for area plotting)
    maxf(i)=sum(f(peakdex).*(nex2(peakdex,i)-thresh)')/sum(nex2(peakdex,i)-thresh); %compute the expected frequency value fo directon i
    centerdex=find(abs(f-maxf(i))==min(abs(f-maxf(i)))); %index of the frequency bin in 'f' nearest to the expected frequency
    %sharpness(i) = sum(nex2((centerdex-8):(centerdex+8),i))/sum(nex2((centerdex-25):(centerdex+25),i)); %measure of how good the theta peak is
end

%% ----- plot the power spectra ------
figure(20); clf; %figure for plotting the power spectra in each direction
for i=1:8
    subplot(8,1,i); hold off;
    area(f(rangedex),[nexchopped(rangedex,i) nextops(rangedex,i)]); hold on; line([maxf(i) maxf(i)], [0 20]); set(gca,'XLim',[4 12]);
end

%%----------------------------------------------------
%%  cosine fitting of the burst frequency tuning curve
%%----------------------------------------------------



    dirfunc=[maxf([1:8]) maxf(1)]; %burst frequency tuning curve
    
    %%compute initial parameter estimates for passing to fitter
    [peak, peakoff]=max(dirfunc);  %initial estimate for size and position of cosine peak
    [valley valloff]=min(dirfunc); %initial estimate for size and position of cosine valley
    figure(1+figbase); %plot burst frequency tuning curve in this figure
    %dirfunc(badpeaks)=NaN;
    [pp, base, phase, gof]=cosfit8(x,dirfunc,(peak-valley),peakoff,valley+(peak-valley)/2,[1 1 1 1 1 1 1 1 0]); %do the cosine fit
    legend('off'); %remove legend from the plot

    pdex=round((2*pi-phase)/.7854); %%index nearest to cosine peak
    if (pdex>8)
        pdex=pdex-8;
    end
    vdex=pdex+4;
    if (vdex>8)
        vdex=vdex-8;
    end

    mdgof=gof;


    theta=rad2deg(phase);
    omega=base;

    goodfit(figbase+1)=gof.rsquare;
    prefd(figbase+1) = 2*pi-phase;


% ----- plot the power spectra ------
figure(20); clf; %figure for plotting the power spectra in each direction
for i=1:8
    subplot(8,1,i); hold off;
    area(f(rangedex),[nexchopped(rangedex,i) nextops(rangedex,i)]); hold on; line([maxf(i) maxf(i)], [0 20]); set(gca,'XLim',[6 10]);
end

figure(21); clf;
subplot(8,1,1); bar(-blocksec:(blocksec/256):blocksec,sum(e_autohisto),1); axis tight;
subplot(8,1,2); bar(-blocksec:(blocksec/256):blocksec,sum(ne_autohisto),1); axis tight;
subplot(8,1,3); bar(-blocksec:(blocksec/256):blocksec,sum(n_autohisto),1); axis tight;
subplot(8,1,4); bar(-blocksec:(blocksec/256):blocksec,sum(nw_autohisto),1); axis tight;
subplot(8,1,5); bar(-blocksec:(blocksec/256):blocksec,sum(w_autohisto),1); axis tight;
subplot(8,1,6); bar(-blocksec:(blocksec/256):blocksec,sum(sw_autohisto),1); axis tight;
subplot(8,1,7); bar(-blocksec:(blocksec/256):blocksec,sum(s_autohisto),1); axis tight;
subplot(8,1,8); bar(-blocksec:(blocksec/256):blocksec,sum(se_autohisto),1); axis tight;
allauto=mean([sum(e_autohisto)/sum(e_ratehisto(:)); sum(ne_autohisto)/sum(ne_ratehisto(:)); sum(n_autohisto)/sum(n_ratehisto(:)); sum(nw_autohisto)/sum(nw_ratehisto(:)); sum(w_autohisto)/sum(w_ratehisto(:)); sum(sw_autohisto)/sum(sw_ratehisto(:)); sum(s_autohisto)/sum(s_ratehisto(:)); sum(se_autohisto)/sum(se_ratehisto(:))]);
figure(19);  bar(-blocksec:(blocksec/256):blocksec,allauto*640,1); axis tight;







