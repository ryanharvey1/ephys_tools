function  [CellClass] = bz_CellClassification_adapted (basepath)
% Loads cell spike waveforms from the local folder, characterizes them and
% separates them into E vs I cells.  Manual verification based on clickable
% gui from Shige.
% 
%INPUTS
% basepath: path to processed data folder
%
%OUTPUTS
% 
% Mixture of functions from Shigeyoshi Fujisawa (WaveShapeClassification),
% Adrien Peyrache (wavelet-based determination of spike width) and Eran Stark (wfeatures, spikestats).
%
% Brendon Watson 2014
% Modified to buzcode format DLevenstein 2017
% Modified to ephys_tools format by L.Berkowitz 2019

%Done(Berkowitz 2/19/20): 
%       - Loads ephys_tools data structure and pulls highest amplitude
%       waveform for classificaiton. 
%       - calulates peak-trough and trough-peak time (in ms) and
%       spike-duration using peyrache wavelet. 
%       - plots prelim classification results. 

%To-do (Berkowitz 2/19/20): 
%       - Add additional features to improve clustering (e.g. Burst idx)
%       - Decide which classification approach to use. K-means is used now,
%       but requires tinkering given the current features used. 
%       - Testing/validation of classificaiton approach. Perhaps use
%       Buzaski lab waveforms that have already been classified? 
%       - Create method of saving waveform classification back to data structure (or separately). 
%       - Create publication quality figure showing classificaiton results
%       :) 
%

%%
sessions = dir([basepath,'*.mat']);
all_waves = [];
wave_id = [];
for s = 1:length(sessions)
    temp = load(fullfile(sessions(s).folder,sessions(s).name),'avgwave','spikesID');
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
    
    thiswave = MaxWaves(a,:); keep_wave(a,:) = thiswave; 
    [~,minpos] = min(thiswave); %find the trough
    minpos = minpos(1);
    
    [~,maxpos] = max(thiswave); %find the largest peak 
    maxpos = maxpos(1); 
    
    if abs(thiswave(maxpos)) > abs(thiswave(minpos))
        thiswave = thiswave*-1; %flip waveform so peak becomes trough - for waveforms that have large trough preceeding peak.
        [~,minpos] = min(thiswave); %recalculate trough
        [~,maxpos] = max(thiswave); %recalculate peak
        minpos = minpos(1);
    end
      
    if maxpos >= 25 || isempty(maxpos) || isempty(minpos) || minpos >= 31 || minpos <= 5
        warning('Your Waveform may be erroneous')
        peak_trough(a,:) = NaN;
        trough_peak(a,:) = NaN;
        continue
    end
    % Find time from min to left
    [max_left,~] = max(thiswave(1,1:minpos-1)); %finds peak before trough
    maxpos_left = find(thiswave == max_left(1));
    
    peak_trough(a,:) = (minpos - maxpos_left(1))/32; % in ms
    
    % Find time from min to right
    [max_right,~] = max(thiswave(1,minpos+1:end)); %finds peak after trough
    maxpos_right = find(thiswave == max_right(1));
    
    trough_peak(a,:) = (maxpos_right(1) -  minpos)/32; % in ms
%     
%     figure; 
%     plot(thiswave); 
%     hold on; 
%     plot(minpos,thiswave(minpos),'*b'); hold on;
%     plot(maxpos,thiswave(maxpos),'*r'); hold on; 
%     plot(maxpos_right,thiswave(maxpos_right),'*k');
%     plot(maxpos_left,thiswave(maxpos_left),'*c');
%     title(['peak2trough =  ', num2str(peak_trough(a,:)),'  trough2peak =  ', num2str(trough_peak(a,:))])
%     close all
end

%% get spike width by taking inverse of max frequency in spectrum (based on Adrien's use of Eran's getWavelet)
for a = 1:size(MaxWaves,1)
    w = MaxWaves(a,:)';
    [~,maxpos] = max(w); %find the largest peak
    maxpos = maxpos(1);
     if maxpos >= 25 || isempty(maxpos) || isempty(minpos) || minpos >= 31
        warning('Your Waveform may be erroneous')
        spkW(a,1) = NaN;
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
plot(mean(keep_wave(idx ==2,:),1),'b')
hold on;
plot(mean(keep_wave(idx ==1,:),1),'r')

x = trough_peak;
y = spkW;%width in ms of wavelet representing largest feature of spike complex... ie the full trough including to the tip of the peak

xx = [0 0.8];
yy = [2.4 0.4];
m = diff( yy ) / diff( xx );
b = yy( 1 ) - m * xx( 1 );  % y = ax+b
RS = y>= m*x+b;
INT = ~RS;


end
