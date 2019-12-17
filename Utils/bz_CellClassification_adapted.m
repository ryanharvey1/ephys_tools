function  [CellClass] = bz_CellClassification_adapted (basepath)
% Loads cell spike waveforms from the local folder, characterizes them and
% separates them into E vs I cells.  Manual verification based on clickable
% gui from Shige.
% 
%INPUTS
% baseName - basename of files in the local folder (default: pwd)
% 'knownE' - UIDs of known E cells, ie from synaptic interactions
% 'knownI' - UIDs of known I cells, ie from synaptic interactions
% 'saveMat'- true/false, save basePath/baseName.CellClass.cellinfo.mat
%            (default:true)
% 'saveFig'- true/false, save a DetectionFigure for posterity/QC 
%            (default:true)
%
%OUTPUTS
%   CellClass   buzcode structure saved to
%               basePath/baseName.CellClass.cellinfo.mat
%       .UID    -UID for each of the cells, matching spikes.cellinfo.mat
%       .pE 	-index vector, true for putative excitatory (RS) cells
%       .pI     -index vector, true for putative inhibitory (NS) cells
%       .label 	-labels for each cell 'pE' or 'pI'
%       .detectionparms.Waveforms -mean waveforms of each cell at the max channel
%       .detectionparms.TroughPeakMs
%       .detectionparms.SpikeWidthMs
%       .detectionparms.PyrBoundary - x,y of manually drawn line of boundary
%
%
% Mixture of functions from Shigeyoshi Fujisawa (WaveShapeClassification),
% Adrien Peyrache (wavelet-based determination of spike width) and Eran Stark (wfeatures, spikestats).
%
% Brendon Watson 2014
% Modified to buzcode format DLevenstein 2017
% Modified to ephys_tools format by L.Berkowitz 2019

%% input Parsing

p = inputParser;
addParameter(p,'knownE',[],@isvector);
addParameter(p,'knownI',[],@isvector);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveFig',true,@islogical);
addParameter(p,'forceReload',false,@islogical);

parse(p,varargin{:})

knownE = p.Results.knownE;
knownI = p.Results.knownI;
SAVEMAT = p.Results.saveMat;



%%
sessions = dir([basepath,'*.mat']);
all_waves = [];
wave_id = [];
for s = 1:length(sessions)
    temp = load(sessions(s).name,'avgwave','spikesID');
    wave_id = [wave_id; cellstr(repmat(sessions(s).name,size(temp.avgwave,1),1))...
        temp.spikesID.TetrodeNum num2cell(temp.spikesID.CellNum)];
    all_waves = [all_waves; temp.avgwave]; 
end

%% gather waves & select waveform with maximum amplitude
for wave = 1:size(all_waves,1) 
    
    %find highest amplitude waveform. 
    temp_waves = all_waves{wave}; 
    [~,maxpos] = max(max(abs(temp_waves ),[],2)); %use abs to capture peak amplitude, positive or negative, whichever is greater.
    wf = temp_waves (maxpos,:);
    
    %because of differences in spike sorting methods, stored waveforms are
    %different lengths. Let's scale to 1 ms or 32 samples. 
    wf_x = linspace(1,32,length(wf));
    wf_32 = interp1(wf_x,wf,1:32);
    
    MaxWaves(wave,:) = wf_32;
    
    clear wf_32 wf_x wf maxpos temp_waves
end

clear wave
%% get trough-peak delay times
for a = 1:size(MaxWaves,1)
    
    thiswave = MaxWaves(a,:);
    [~,minpos] = min(thiswave); %find the trough
    minpos = minpos(1);
    
    [~,maxpos] = max(thiswave); %find the largest peak 
    maxpos = maxpos(1); 
    
    if maxpos >= 25 || isempty(maxpos) || isempty(minpos) || minpos >= 31
        warning('Your Waveform may be erroneous')
        peak_trough(a,:) = 4;
        trough_peak(a,:) = 4;
        continue
    end
    
    if minpos < maxpos
        thiswave = thiswave*-1; %flip waveform so peak becomes trough - for waveforms that have large trough preceeding peak.
        [~,minpos] = min(thiswave); %recalculate trough
        minpos = minpos(1);
    end
    
    % Find time from min to left
    [max_left,~] = max(thiswave(1,1:minpos-1)); %finds peak before trough
    maxpos_left = find(thiswave == max_left(1));
    
    peak_trough(a,:) = (minpos - maxpos_left(1))/32; % in ms
    
    % Find time from min to right
    [max_right,~] = max(thiswave(1,minpos+1:end)); %finds peak after trough
    maxpos_right = find(thiswave == max_right(1));
    
    trough_peak(a,:) = (maxpos_right(1) -  minpos)/32; % in ms
    
end

%% get spike width by taking inverse of max frequency in spectrum (based on Adrien's use of Eran's getWavelet)
for a = 1:size(MaxWaves,1)
    w = MaxWaves(a,:)';
    [~,maxpos] = max(w); %find the largest peak
    maxpos = maxpos(1);
     if maxpos >= 25 || isempty(maxpos) || isempty(minpos) || minpos >= 31
        warning('Your Waveform may be erroneous')
        spkW(a,1) = 4;
        continue
    end
    w = [w(1)*ones(1000,1);w;w(end)*ones(1000,1)];
    [wave, f, t] = getWavelet(w,32000,300,8000,128);
    %We consider only the central portion of the wavelet because we
    %haven't filtered it before hand (e.g. with a Hanning window)
    wave = wave(:,int16(length(t)/4):3*int16(length(t)/4));
    %Where is the max frequency?
    [maxPow, ix] = max(wave);
    [~, mix] = max(maxPow);
    ix = ix(mix);
    spkW(a,1) = 1000/f(ix);
end

clear wave 

%% Generate separatrix for cells 
X=[trough_peak,spkW];
[idx,C] = kmeans(X,2, 'Replicates',10,'Distance','cityblock');

figure;
plot(X(idx==1,1),X(idx==1,2),'r.','MarkerSize',12)
hold on
plot(X(idx==2,1),X(idx==2,2),'b.','MarkerSize',12)
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
title 'Cluster Assignments and Centroids'
hold off

figure; 
plot(mean(MaxWaves(idx ==2,:),1),'r')
hold on;
plot(mean(MaxWaves(idx ==1,:),1),'b')

x = trough_peak;
y = spkW;%width in ms of wavelet representing largest feature of spike complex... ie the full trough including to the tip of the peak

xx = [0 0.8];
yy = [2.4 0.4];
m = diff( yy ) / diff( xx );
b = yy( 1 ) - m * xx( 1 );  % y = ax+b
RS = y>= m*x+b;
INT = ~RS;


end
