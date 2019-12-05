function  [CellClass] = bz_CellClassification (data)
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
OneMs = round(32000/1000); %sample rate is 32k hz. divide by 1000 to put into ms

%% gather waves & select waveform with maximum amplitude
for wave = 1:size(data.avgwave,1) 
    
    %find highest amplitude waveform. 
    all_waves = data.avgwave{wave}; 
    [~,maxpos] = max(max(abs(all_waves),[],2)); %use abs to capture peak amplitude, positive or negative, whichever is greater.
    wf = all_waves(maxpos,:);
    
    %because of differences in spike sorting methods, stored waveforms are
    %different lengths. Let's scale to 1 ms or 32 samples. 
    wf_x = linspace(1,32,length(wf));
    wf_32 = interp1(wf_x,wf,1:32);
    
    MaxWaves(wave,:) = wf_32;
end

%% get trough-peak delay times
for a = 1:size(MaxWaves,1)
    thiswave = MaxWaves(a,:);
    
    
    [minval,minpos] = min(thiswave);
    minpos = minpos(1);
    [maxval,maxpos] = max(thiswave);
        [dummy,maxpos] = max(thiswave(minpos+1:end));
        if isempty(maxpos)
            warning('Your Waveform may be erroneous')
            maxpos = 1
        end
        maxpos=maxpos(1);
        maxpos = maxpos+minpos;
        tp(a) = maxpos-minpos; %In number of samples
end

%% get spike width by taking inverse of max frequency in spectrum (based on Adrien's use of Eran's getWavelet)
for a = 1:size(MaxWaves,2)
    w = MaxWaves(:,a);
    w = [w(1)*ones(1000,1);w;w(end)*ones(1000,1)];
    [wave f t] = getWavelet(w,20000,500,3000,128);
    %We consider only the central portion of the wavelet because we
    %haven't filtered it before hand (e.g. with a Hanning window)
    wave = wave(:,int16(length(t)/4):3*int16(length(t)/4));
    %Where is the max frequency?
    [maxPow ix] = max(wave);
    [dumy mix] = max(maxPow);
    ix = ix(mix);
    spkW(a) = 1000/f(ix);
end

%% Generate separatrix for cells 
x = tp'/OneMs;%trough to peak in ms
y = spkW';%width in ms of wavelet representing largest feature of spike complex... ie the full trough including to the tip of the peak

xx = [0 0.8];
yy = [2.4 0.4];
m = diff( yy ) / diff( xx );
b = yy( 1 ) - m * xx( 1 );  % y = ax+b
RS = y>= m*x+b;
INT = ~RS;

%% Plot for manual selection of boundary, with display of separatrix as a guide.
h = figure;
title({'Discriminate pyr and int (select Pyramidal)','left click to draw boundary', 'center click/ENTER to complete)'});
fprintf('\nDiscriminate pyr and int (select Pyramidal)');
xlabel('Trough-To-Peak Time (ms)')
ylabel('Wave width (via inverse frequency) (ms)')
[ELike,PyrBoundary] = ClusterPointsBoundaryOutBW([x y],knownE,knownI,m,b);

%% Mean waveforms output
CellClass.UID = spikes.UID; %ADD CELL IDX 
CellClass.pE = ELike';
CellClass.pI = ~ELike';
CellClass.label = cell(size(CellClass.UID));
CellClass.label(CellClass.pE) = {'pE'};
CellClass.label(CellClass.pI) = {'pI'};
CellClass.detectionparms.TroughPeakMs = x';
CellClass.detectionparms.SpikeWidthMs = y';
CellClass.detectionparms.PyrBoundary = PyrBoundary;
CellClass.detectionparms.Waveforms = MaxWaves;

if SAVEMAT
    save(savefile,'CellClass')
end

%%
if SAVEFIG
    figure
    subplot(2,2,1)
        plot(CellClass.detectionparms.TroughPeakMs(CellClass.pE),...
            CellClass.detectionparms.SpikeWidthMs(CellClass.pE),'k.')
        hold on
        plot(CellClass.detectionparms.TroughPeakMs(CellClass.pI),...
            CellClass.detectionparms.SpikeWidthMs(CellClass.pI),'r.')
        axis tight
        plot(CellClass.detectionparms.PyrBoundary(:,1),...
            CellClass.detectionparms.PyrBoundary(:,2))
        xlim([0 max([x+0.1;2])])
        ylim([0 max([y+0.1;2])])
        xb = get(gca,'XLim');
        yb = get(gca,'YLim');
        plot(xb,[m*xb(1)+b m*xb(2)+b])
        xlabel('Trough to Peak Time (ms)')
        ylabel('Spike Width (ms)')
        
    subplot(2,2,2)
        plot([1:size(MaxWaves,1)]./OneMs,MaxWaves(:,CellClass.pE),'color',[0 0.6 0])
        hold on
        plot([1:size(MaxWaves,1)]./OneMs,MaxWaves(:,CellClass.pI),'color',[0.6 0 0])
        axis tight
        xlabel('t (ms)')
        
end
