function  [CellClass] = bz_CellClassification_adapted(basepath,varargin)
% Loads cell spike waveforms from the local folder, characterizes them and
% separates them into E vs I cells.  Manual verification based on clickable
% gui from Shige.
%
%INPUTS
% basepath: path to project folder
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
%Done(Harvey 2/20/20)
%       - Add additional features to improve clustering
%       - Create method of saving waveform classification back to data structure

%To-do (Berkowitz 2/19/20):
%       - Decide which classification approach to use. K-means is used now,
%       but requires tinkering given the current features used.
%       - Testing/validation of classificaiton approach. Perhaps use
%       Buzaski lab waveforms that have already been classified?
%       - Create publication quality figure showing classificaiton results
%       :)
%
% dependencies:
%   getWavelet: ~\ephys_tools\external_packages\buzcode\analysis\spikes\cellTypeClassification\BrendonClassificationFromStark2013
%   WaveletToolbox: ~\ephys_tools\external_packages\buzcode\externalPackages\WaveletToolbox
%   CCG: ~\ephys_tools\external_packages\buzcode\analysis\spikes\correlation
%   cell-explorer https://github.com/petersenpeter/Cell-Explorer
%%
p = inputParser;
addParameter(p,'manualAdjustMonoSyn',0);
addParameter(p,'load_saved_cell_metrics',0);
parse(p,varargin{:})
manualAdjustMonoSyn = p.Results.manualAdjustMonoSyn;
load_saved_cell_metrics = p.Results.load_saved_cell_metrics;


if load_saved_cell_metrics
    %check for cell metrics before compiling features
    if exist(fullfile(basepath,'cell_metrics.mat'),'file')
        load(fullfile(basepath,'cell_metrics.mat'),'cell_metrics')
    end
end

if ~exist('cell_metrics','var')
    
    sessions = dir([basepath,'\ProcessedData\*.mat']);
    all_waves = {};
    wave_id = {};
    cell_metrics = [];
    
    % get n cells per session
    i = 0;
    for s = 1:length(sessions)
        load(fullfile(sessions(s).folder,sessions(s).name),'spikesID');
        idx{s} = i+1:length(spikesID.TetrodeNum)+i;
        i = i+length(spikesID.TetrodeNum);
    end
    
    WaitMessage = parfor_wait(length(sessions));
    
    for s = 1:length(sessions)
        temp = load(fullfile(sessions(s).folder,sessions(s).name),...
            'avgwave','spikesID','Spikes','session_path','rat','sessionID');
        
        % autocor metrics
        [acg_metrics] = calc_acg_metrics(temp);
        fit_params = fit_ACG(acg_metrics.acg2);
        
        % compile metrics
        cell_metrics = get_cell_metrics(cell_metrics,acg_metrics,fit_params,idx,s);
        
        % get cell id
        cell_metrics.wave_id(idx{s},:) = [cellstr(repmat(sessions(s).name,size(temp.avgwave,1),1)),...
            temp.spikesID.TetrodeNum, num2cell(temp.spikesID.CellNum)];
        
        % standardize into 32 samples & pull out largest amp channel
        cell_metrics.waveforms(idx{s},:) = get_waveform(temp);
        
        % Pytative MonoSynaptic connections
        [cell_metrics] = mono_synaptic(temp,cell_metrics,s,manualAdjustMonoSyn); 
        
        WaitMessage.Send;
    end
    WaitMessage.Destroy;
    
    % get trough-peak delay times and flip waveforms that need flipped
    [cell_metrics] = get_waveform_delay_times(cell_metrics);
    
    % get spike width by taking inverse of max frequency in spectrum
    % (based on Adrien's use of Eran's getWavelet)
    cell_metrics = get_spike_wavelet(cell_metrics);
    
    % zscore waveforms
    cell_metrics.waveforms_zscore = zscore(cell_metrics.waveforms,0,2);
    
    % align waveforms (trough at the 8th column)
    cell_metrics = align_waveforms(cell_metrics);

    % save features
    save(fullfile(basepath,'cell_metrics.mat'),'cell_metrics')
end

%% Generate separatrix for cells
 X = [cell_metrics.thetaModulationIndex,...
    cell_metrics.burstIndex_Royer2012,...
    cell_metrics.burstIndex_Doublets,...
    cell_metrics.acg_tau_decay,...
    cell_metrics.acg_tau_rise,...
    cell_metrics.acg_c,...
    cell_metrics.acg_d,...
    cell_metrics.acg_asymptote,...
    cell_metrics.acg_refrac,...
    cell_metrics.acg_fit_rsquare,...
    cell_metrics.acg_tau_burst,...
    cell_metrics.acg_h,...
    cell_metrics.peakToTrough,...
    cell_metrics.troughToPeak,...
    cell_metrics.spkW];
X(isinf(X)) = NaN;

X = normalize(X,1,'range');
% X = normalize(X,1);
% [pc,score,latent,tsquare] = pca(X);
% X = score(:,1);
[idx,C] = kmeans(X,2, 'Replicates',10,'Distance','cityblock');


figure;
subplot(2,2,1)
scatter(cell_metrics.troughToPeak(idx==1),cell_metrics.spkW(idx==1),...
    12,'r','Filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
hold on
scatter(cell_metrics.troughToPeak(idx==2),cell_metrics.spkW(idx==2),...
    12,'b','Filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
% plot(C(:,1),C(:,2),'kx',...
%     'MarkerSize',15,'LineWidth',3)
legend('Cluster 1','Cluster 2','Centroids',...
    'Location','NW')
title 'Cluster Assignments and Centroids'
hold off

% classify cell type by shortest trough to peak time
[~,I] = min([nanmean(cell_metrics.troughToPeak(idx == 1)),...
    nanmean(cell_metrics.troughToPeak(idx == 2))]);
cell_type = {'Interneuron','Pyr'};


subplot(2,2,3)
plot(1:32,cell_metrics.aligned_waves(idx == 2,:),'Color',[0,0,1,0.1]);
hold on;
plot(nanmean(cell_metrics.aligned_waves(idx == 2,:),1),'k','Linewidth',3)
ylabel('zscore')
xlabel('ms')
if I == 2
    title(cell_type{1})
else
    title(cell_type{2})
end
subplot(2,2,4)
plot(1:32,cell_metrics.aligned_waves(idx == 1,:),'Color',[1,0,0,0.1]);
hold on;
plot(nanmean(cell_metrics.aligned_waves(idx == 1,:),1),'k','Linewidth',3)
if I == 1
    title(cell_type{1})
else
    title(cell_type{2})
end



% below 0.25 ms
idx = cell_metrics.troughToPeak <= 0.25;
figure;
subplot(2,2,1)
scatter(X(idx,1),X(idx,2),12,'r','Filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
hold on
scatter(X(~idx,1),X(~idx,2),12,'b','Filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

subplot(2,2,3)
p1 = plot(1:32,cell_metrics.aligned_waves(idx,:),'Color',[0,0,1,0.1]);
hold on;
plot(nanmean(cell_metrics.aligned_waves(idx,:),1),'k','Linewidth',3)
ylabel('zscore')
xlabel('ms')
title(cell_type{1})
subplot(2,2,4)
p2 = plot(1:32,cell_metrics.aligned_waves(~idx,:),'Color',[1,0,0,0.1]);
hold on;
plot(nanmean(cell_metrics.aligned_waves(~idx,:),1),'k','Linewidth',3)
title(cell_type{2})



cell_metrics.filled_waves = fillmissing(cell_metrics.aligned_waves,'previous',2,'EndValues','nearest');
[pc,score,latent,tsquare] = pca(cell_metrics.filled_waves);

eva = evalclusters(score(:,1:3),'kmeans','CalinskiHarabasz','KList',[1:6]);

% K-mean clustering
opts = statset('Display','final');
klusters = eva.OptimalK;
colors = {'r','b','g','k'};
colors2 = [1,0,0,0.1;0,1,0,0.1;0,0,1,0.1;0,0,1,0.1]';
[idx,C] = kmeans(score(:,1:3),klusters,'Distance','cityblock',...
    'Replicates',5,'Options',opts);
% waveform_metrics.klusters = idx';

figure
subplot(2,2,1:2)
for i = 1:klusters
    plot(cell_metrics.filled_waves(idx==i,:)','Color',colors2(:,i),'linewidth',2), hold on
end
xlabel('Time (ms)'),ylabel('Z-scored amplitude'),title 'Average Waveforms', axis tight

subplot(2,2,3)
for i = 1:klusters
    scatter3(score(idx==i,1)',score(idx==i,2)',score(idx==i,3)',20,colors{i}), hold on
end
plot3(C(:,1),C(:,2),C(:,3),'kx','MarkerSize',15,'LineWidth',3), title 'PCA Analysis', hold off
xlabel('PC1'),ylabel('PC2'),zlabel('PC3'),axis tight

for i = 1:klusters
    figure
    plot(cell_metrics.filled_waves(idx==i,:)','Color',colors2(:,i),'linewidth',2);
    xlabel('Time (ms)'),ylabel('Z-scored amplitude'),title 'Average Waveforms', axis tight
end




% cell_classification_putativeCellType
cell_metrics.general.cellCount = length(cell_metrics.thetaModulationIndex);


%     dispLog('Performing Cell-type classification');
% All cells assigned as Pyramidal cells at first
cell_metrics.putativeCellType = repmat({'Pyramidal Cell'},1,cell_metrics.general.cellCount);

% Interneuron classification
% Cells are reassigned as interneurons by below criteria
% acg_tau_decay > 30ms
cell_metrics.putativeCellType(cell_metrics.acg_tau_decay>30) = ...
    repmat({'Interneuron'},sum(cell_metrics.acg_tau_decay>30),1);
% acg_tau_rise > 3ms
cell_metrics.putativeCellType(cell_metrics.acg_tau_rise>3) = ...
    repmat({'Interneuron'},sum(cell_metrics.acg_tau_rise>3),1);
% troughToPeak <= 0.25ms
cell_metrics.putativeCellType(cell_metrics.troughToPeak<=0.25) = ...
    repmat({'Interneuron'},sum(cell_metrics.troughToPeak<=0.25),1);
% Narrow interneuron assigned if troughToPeak <= 0.25ms
cell_metrics.putativeCellType([cell_metrics.troughToPeak<=0.25]' & ...
    ismember(cell_metrics.putativeCellType, 'Interneuron')) = ...
    repmat({'Narrow Interneuron'},sum([cell_metrics.troughToPeak<=0.25]' &...
    (ismember(cell_metrics.putativeCellType, 'Interneuron'))),1);
% Wide interneuron assigned if troughToPeak > 0.25ms
cell_metrics.putativeCellType([cell_metrics.troughToPeak>0.25]' &...
    ismember(cell_metrics.putativeCellType, 'Interneuron')) = ...
    repmat({'Wide Interneuron'},sum([cell_metrics.troughToPeak>0.25]' &...
    (ismember(cell_metrics.putativeCellType, 'Interneuron'))),1);

figure
colors = {'r','b','g'};
colors2 = [1,0,0,0.1;0,1,0,0.1;0,0,1,0.1]';

subplot(1,4,1)
plot(cell_metrics.filled_waves(cell_metrics.putativeCellType == "Pyramidal Cell",:)','Color',colors2(:,1),'linewidth',2);
xlabel('Time (ms)'),ylabel('Z-scored amplitude'),title 'Average Waveforms', axis tight
hold on
plot(cell_metrics.filled_waves(cell_metrics.putativeCellType == "Wide Interneuron",:)','Color',colors2(:,2),'linewidth',2);
xlabel('Time (ms)'),ylabel('Z-scored amplitude'),title 'Average Waveforms', axis tight

plot(cell_metrics.filled_waves(cell_metrics.putativeCellType == "Narrow Interneuron",:)','Color',colors2(:,3),'linewidth',2);
xlabel('Time (ms)'),ylabel('Z-scored amplitude'),title 'Average Waveforms', axis tight


subplot(1,4,2)
plot(cell_metrics.filled_waves(cell_metrics.putativeCellType == "Pyramidal Cell",:)','Color',colors2(:,1),'linewidth',2);
xlabel('Time (ms)'),ylabel('Z-scored amplitude'),title 'Average Waveforms', axis tight
title('Pyramidal Cell')

subplot(1,4,3)
plot(cell_metrics.filled_waves(cell_metrics.putativeCellType == "Wide Interneuron",:)','Color',colors2(:,2),'linewidth',2);
xlabel('Time (ms)'),ylabel('Z-scored amplitude'),title 'Average Waveforms', axis tight
title('Wide Interneuron')

subplot(1,4,4)
plot(cell_metrics.filled_waves(cell_metrics.putativeCellType == "Narrow Interneuron",:)','Color',colors2(:,3),'linewidth',2);
xlabel('Time (ms)'),ylabel('Z-scored amplitude'),title 'Average Waveforms', axis tight
title('Narrow Interneuron')

figure;

subplot(1,3,1)
imagesc(zscore(cell_metrics.acg.narrow(:,cell_metrics.putativeCellType == "Pyramidal Cell")',[],2))
title('Pyramidal Cell')
ax = gca;
mag_range(1,:) = ax.CLim;
set(ax,'XTick',linspace(1,201,5),'XTickLabel',linspace(-100,100,5))
xlabel('lag (ms)')
ylabel('Unit')

subplot(1,3,2)
imagesc(zscore(cell_metrics.acg.narrow(:,cell_metrics.putativeCellType == "Wide Interneuron")',[],2))
title('Wide Interneuron')
ax = gca;
mag_range(2,:) = ax.CLim;
set(ax,'XTickLabel',[])

subplot(1,3,3)
imagesc(zscore(cell_metrics.acg.narrow(:,cell_metrics.putativeCellType == "Narrow Interneuron")',[],2))
title('Narrow Interneuron')
ax = gca;
mag_range(3,:) = ax.CLim;
set(ax,'XTickLabel',[])


sgtitle('Narrow acg') 

for i = 1:3
subplot(1,3,i)
ax = gca;
ax.CLim = [min(mag_range(:)), max(mag_range(:))];
end
colormap(viridis(255))

end

function [acg_metrics] = calc_acg_metrics(temp)

[ccg,time] = CCG(temp.Spikes,[],'binSize',0.001,'duration',1,'norm','rate','Fs',32000);
[ccg2,time2] = CCG(temp.Spikes,[],'binSize',0.0005,'duration',0.100,'norm','rate','Fs',32000);
max_lags = (length(time)-1)/2;
max_lags2 = (length(time2)-1)/2;
for ii = 1:size(ccg,2)
    xc = ccg(:,ii,ii);
    xc2 = ccg2(:,ii,ii);
    ThetaModulationIndex(ii) = (mean(xc(max_lags+1+50:max_lags+1+70))-...
        mean(xc(max_lags+1+100:max_lags+1+140)))/...
        (mean(xc(max_lags+1+50:max_lags+1+70))+...
        mean(xc(max_lags+1+100:max_lags+1+140)));
    BurstIndex_Royer2012(ii,1) = sum(xc(max_lags+1+3:max_lags+1+5))/...
        mean(xc(max_lags+1+200:max_lags+1+300));
    BurstIndex_Doublets(ii,1) = max(xc2(max_lags2+1+5:max_lags2+1+16))/...
        mean(xc2(max_lags2+1+16:max_lags2+1+23));
    ACG(:,ii) = xc;
    ACG2(:,ii) = ccg2(:,ii,ii);
end
acg_metrics.acg = ACG;
acg_metrics.acg2 = ACG2;
acg_metrics.thetaModulationIndex = ThetaModulationIndex;
acg_metrics.burstIndex_Royer2012 = BurstIndex_Royer2012;
acg_metrics.burstIndex_Doublets = BurstIndex_Doublets;
acg_metrics.ccg = ccg2;
acg_metrics.ccg_time = time2;

end

function fit_params_out = fit_ACG(acg2)
figs = false;

fit_params = [];
rsquare = [];
offset = 101;

g = fittype('max(c*(exp(-(x-f)/a)-d*exp(-(x-f)/b))+h*exp(-(x-f)/g)+e,0)',...
    'dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c','d','e','f','g','h'});

acg_mat =[]; paut=[];
jj = 1;
for j = 1:size(acg2,2)
    x = ([1:100]/2)';
    y = acg2(x*2+offset,j); % /max(acg2(x+offset,j))
    
    [f0,gof] = fit(x,y,g,'StartPoint',[20, 1, 30, 2, 0, 5, 1.5,2],...
        'Lower',[1, 0.1, 0, 0, -30, 0,0.1,0],...
        'Upper',[500, 50, 500, 15, 50, 20,5,100]);
    %     f0 = fit(x,y,'exp2','StartPoint',[1, -0.015, -1, -1]);
    fit_params(:,j) = coeffvalues(f0);
    %     xx = linspace(0,48,100);
    rsquare(j) = gof.rsquare;
    if figs
        if rem(j,12)==1
            figure('position',[50,50,1000,800]),
            jj = 1;
        end
        subplot(3,4,jj)
        bar_from_patch(x,y,'b'), hold on
        plot(x,f0(x),'r-');
        axis tight, title([num2str(j) , ': rise=' num2str(fit_params(2,j),3),...
            ', decay=' num2str(fit_params(1,j),3),])
        %     ylim([0,1])
        %     subplot(2,1,2)
        %     [fmodel,ydata,xdata,paut(j,:)] = fitpyrint(acg2(:,j)',0:0.5:50,1,20);
        
        jj = jj + 1;
    end
end
if figs
    figure,
    subplot(3,3,1), x1 = fit_params(1,:);
    [~,edges] = histcounts(log10(x1),40);
    histogram(x1,10.^edges), xlabel('\tau_{decay} [ms]'), axis tight, ylabel('A')
    set(gca, 'xscale','log')
    
    subplot(3,3,2), x1 = fit_params(2,:);
    [~,edges] = histcounts(log10(x1),40);
    histogram(x1,10.^edges), xlabel('\tau_{rise} [ms]'), axis tight,...
        title('double-exponential fit to the ACG'), ylabel('B')
    set(gca, 'xscale','log')
    
    subplot(3,3,3), x1 = fit_params(3,:);
    [~,edges] = histcounts(log10(x1),40);
    histogram(x1,10.^edges), xlabel('constant c'), axis tight, ylabel('C')
    set(gca, 'xscale','log')
    
    subplot(3,3,4), x1 = fit_params(4,:);
    [~,edges] = histcounts(log10(x1),40);
    histogram(x1,10.^edges), xlabel('constant d'), axis tight, ylabel('D')
    set(gca, 'xscale','log')
    
    subplot(3,3,5), x1 = fit_params(5,:);
    histogram(x1,40), xlabel('constant e'), axis tight, ylabel('E')
    
    subplot(3,3,6), x1 = fit_params(6,:);
    histogram(x1,40), xlabel('t_0 Refractory time (ms)'), axis tight, ylabel('F')
    
    subplot(3,3,7), x1 = fit_params(7,:);
    histogram(x1,40), xlabel('\tau_{burst}'), axis tight, ylabel('G')
    
    subplot(3,3,8), x1 = fit_params(8,:);
    histogram(x1,40), xlabel('Burst constant'), axis tight,...
        title('fit = max(c*(exp(-(x-t_0)/\tau_{decay})-d*exp(-(x-t_0)/\tau_{rise}))+h*exp(-(x-t_0)/\tau_{burst})+e,0)'), ylabel('H')
    
    subplot(3,3,9), x1 = rsquare;
    histogram(x1,40), xlabel('r^2'), axis tight, ylabel('r^2')
end
fit_params_out.acg_tau_decay = fit_params(1,:);
fit_params_out.acg_tau_rise = fit_params(2,:);
fit_params_out.acg_c = fit_params(3,:);
fit_params_out.acg_d = fit_params(4,:);
fit_params_out.acg_asymptote = fit_params(5,:);
fit_params_out.acg_refrac = fit_params(6,:);
fit_params_out.acg_tau_burst = fit_params(7,:);
fit_params_out.acg_h = fit_params(8,:);
fit_params_out.acg_fit_rsquare = rsquare;

end

function bar_from_patch(x_data, y_data,col)
x_data = [x_data(1),reshape([x_data,x_data([2:end,end])]',1,[]),x_data(end)];
y_data = [0,reshape([y_data,y_data]',1,[]),0];
patch(x_data, y_data,col,'EdgeColor','none','FaceAlpha',.8,'HitTest','off')
end

function cell_metrics = get_cell_metrics(cell_metrics,acg_metrics,fit_params,idx,s)

cell_metrics.acg.wide(:,idx{s}) = acg_metrics.acg; % Wide: 1000ms wide CCG with 1ms bins
cell_metrics.acg.narrow(:,idx{s}) = acg_metrics.acg2; % Narrow: 100ms wide CCG with 0.5ms bins
cell_metrics.thetaModulationIndex(idx{s},1) = acg_metrics.thetaModulationIndex; % cell_tmi
cell_metrics.burstIndex_Royer2012(idx{s},1) = acg_metrics.burstIndex_Royer2012; % cell_burstRoyer2012
cell_metrics.burstIndex_Doublets(idx{s},1) = acg_metrics.burstIndex_Doublets;

cell_metrics.acg_tau_decay(idx{s},1) = fit_params.acg_tau_decay;
cell_metrics.acg_tau_rise(idx{s},1) = fit_params.acg_tau_rise;
cell_metrics.acg_c(idx{s},1) = fit_params.acg_c;
cell_metrics.acg_d(idx{s},1) = fit_params.acg_d;
cell_metrics.acg_asymptote(idx{s},1) = fit_params.acg_asymptote;
cell_metrics.acg_refrac(idx{s},1) = fit_params.acg_refrac;
cell_metrics.acg_fit_rsquare(idx{s},1) = fit_params.acg_fit_rsquare;
cell_metrics.acg_tau_burst(idx{s},1) = fit_params.acg_tau_burst;
cell_metrics.acg_h(idx{s},1) = fit_params.acg_h;

cell_metrics.general.ccg{s} = acg_metrics.ccg;
cell_metrics.general.ccg_time{s} = acg_metrics.ccg_time;
end

function waves = get_waveform(temp)

for w = 1:length(temp.avgwave)
    temp.avgwave{w} = interp1(linspace(1,32,length(temp.avgwave{w})),vertcat(temp.avgwave{w})',1:32)';
    [~,I] = max(max(abs(temp.avgwave{w}),[],2));
    waves(w,:) = temp.avgwave{w}(I,:);
end
end

function [cell_metrics]=get_waveform_delay_times(cell_metrics)

for a = 1:size(cell_metrics.waveforms,1)
    
    thiswave = cell_metrics.waveforms(a,:);
    
    [~,minpos] = min(thiswave); %find the trough
    minpos = minpos(1);
    
    [~,maxpos] = max(thiswave); %find the largest peak
    maxpos = maxpos(1);
    
    if abs(thiswave(maxpos)) > abs(thiswave(minpos))
        thiswave = thiswave*-1; %flip waveform so peak becomes trough - for waveforms that have large trough preceeding peak.
        [~,minpos] = min(thiswave); %recalculate trough
        minpos = minpos(1);
        [~,maxpos] = max(thiswave); %recalculate peak
        maxpos = maxpos(1);
        
        % save flipped waveform
        cell_metrics.waveforms(a,:) = thiswave;
        
        % save flipped indicator
        cell_metrics.flipped(a,1) = true;
    else
        % save flipped indicator
        cell_metrics.flipped(a,1) = false;
    end
    
    
    if minpos >= 25 || isempty(maxpos) || isempty(minpos) || minpos ==1
        warning('Your Waveform may be erroneous')
        cell_metrics.peakToTrough(a,:) = NaN;
        cell_metrics.troughToPeak(a,:) = NaN;
        
        %         figure;
        %         plot(thiswave);
        
        continue
    end
    % Find time from min to left
    [max_left,~] = max(thiswave(1,1:minpos-1)); %finds peak before trough
    maxpos_left = find(thiswave == max_left(1));
    
    cell_metrics.peakToTrough(a,:) = (minpos - maxpos_left(1))/32; % in ms
    
    % Find time from min to right
    [max_right,~] = max(thiswave(1,minpos+1:end)); %finds peak after trough
    maxpos_right = find(thiswave == max_right(1));
    
    cell_metrics.troughToPeak(a,:) = (maxpos_right(1) -  minpos)/32; % in ms
    
    
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
end

function cell_metrics = get_spike_wavelet(cell_metrics)

for a = 1:size(cell_metrics.waveforms,1)
    w = cell_metrics.waveforms(a,:)';
    if isnan(cell_metrics.peakToTrough(a))
        cell_metrics.spkW(a,1) = NaN;
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
    cell_metrics.spkW(a,1) = 1000/f(ix);
end
end

function cell_metrics = align_waveforms(cell_metrics)
    packed_waves = [NaN(size(cell_metrics.waveforms_zscore,1),32),cell_metrics.waveforms_zscore,...
        NaN(size(cell_metrics.waveforms_zscore,1),32)];
    [~,trough_idx] = nanmin(packed_waves,[],2);
    for i = 1:size(cell_metrics.waveforms_zscore,1)
        cell_metrics.aligned_waves(i,:) = packed_waves(i,trough_idx(i)-7:trough_idx(i)+24);
    end
end

function [cell_metrics] = mono_synaptic(data,cell_metrics,s,manualAdjustMonoSyn) 
% mono_synaptic: get monosynaptic connections based on spike times
% ephys_tools wrapper for ce_MonoSynConvClick bz_PlotMonoSyn
% dependencies: Cell-Explorer (https://github.com/petersenpeter/Cell-Explorer)
%
% Ryan H 2020

cell_metrics.general.cellCount{s} = length(data.Spikes);

[spikestimes,spikeIDs] = get_ce_MonoSynConvClick_input(data);

% set up saving directory
processedpath=strsplit(data.session_path,filesep);
processedpath(end-2:end)=[];
save_dir = fullfile(strjoin(processedpath,filesep),'mono_res');

% check for saved mono results
if exist(fullfile(save_dir,[data.rat,'_',data.sessionID,'.mono_res.cellinfo.mat']),'file')
    disp('  Loading previous detected MonoSynaptic connections')
    load(fullfile(save_dir,[data.rat,'_',data.sessionID,'.mono_res.cellinfo.mat']),'mono_res');
else
    mono_res = ce_MonoSynConvClick(spikeIDs,spikestimes);
end
spikes.total = mono_res.n;

% manually adjust results
if manualAdjustMonoSyn
    mono_res = gui_MonoSyn(mono_res);
end
if ~exist(save_dir,'dir')
    mkdir(save_dir)
end
save(fullfile(save_dir,[data.rat,'_',data.sessionID,'.mono_res.cellinfo.mat']),...
    'mono_res','-v7.3','-nocompression')

if ~isempty(mono_res.sig_con)
    
    cell_metrics.putativeConnections.excitatory{s} = mono_res.sig_con; % Vectors with cell pairs
    
    cell_metrics.putativeConnections.inhibitory{s} = [];
    
    cell_metrics.synapticEffect{s} = repmat({'Unknown'},1,cell_metrics.general.cellCount{s});
    
    % cell_synapticeffect ['Inhibitory','Excitatory','Unknown']
    cell_metrics.synapticEffect{s}(cell_metrics.putativeConnections.excitatory{s}(:,1)) =...
        repmat({'Excitatory'},1,size(cell_metrics.putativeConnections.excitatory{s},1)); 
    
    cell_metrics.synapticConnectionsOut{s} = zeros(1,cell_metrics.general.cellCount{s});
    
    cell_metrics.synapticConnectionsIn{s} = zeros(1,cell_metrics.general.cellCount{s});
    
    [a,b]=hist(cell_metrics.putativeConnections.excitatory{s}(:,1),...
        unique(cell_metrics.putativeConnections.excitatory{s}(:,1)));
    
    cell_metrics.synapticConnectionsOut{s}(b) = a;
    
    cell_metrics.synapticConnectionsOut{s} = ...
        cell_metrics.synapticConnectionsOut{s}(1:cell_metrics.general.cellCount{s});
    
    [a,b]=hist(cell_metrics.putativeConnections.excitatory{s}(:,2),...
        unique(cell_metrics.putativeConnections.excitatory{s}(:,2)));
    
    cell_metrics.synapticConnectionsIn{s}(b) = a;
    
    cell_metrics.synapticConnectionsIn{s} = ...
        cell_metrics.synapticConnectionsIn{s}(1:cell_metrics.general.cellCount{s});
    
    % Connection strength
    disp('  Determining MonoSynaptic connection strengths (transmission probabilities)')
    for i = 1:size(mono_res.sig_con,1)
        rawCCG = round(cell_metrics.general.ccg{s}(:,mono_res.sig_con(i,1),...
            mono_res.sig_con(i,2))*spikes.total(mono_res.sig_con(i,1))*0.001);
        
        [trans,prob,prob_uncor,pred] = ...
            ce_GetTransProb(rawCCG,spikes.total(mono_res.sig_con(i,1)),0.001,0.020);
        
        cell_metrics.putativeConnections.excitatoryTransProb{s}(i) = trans;
    end
else
    cell_metrics.putativeConnections.excitatory{s} = [];
    cell_metrics.putativeConnections.inhibitory{s} = [];
    cell_metrics.synapticConnectionsOut{s} = zeros(1,cell_metrics.general.cellCount{s});
    cell_metrics.synapticConnectionsIn{s} = zeros(1,cell_metrics.general.cellCount{s});
end
end

function[spikestimes,spikeIDs] = get_ce_MonoSynConvClick_input(data)
spikestimes = vertcat(data.Spikes{:});
i = 0;
for s = 1:length(data.Spikes)
    idx = i+1:length(data.Spikes{s})+i;
    i = i+length(data.Spikes{s});
    spikeIDs(idx,1) =...
        repmat(str2double(extractBetween(data.spikesID.TetrodeNum(s),'TT','.mat')),...
        length(data.Spikes{s}),1);
    spikeIDs(idx,2) = repmat(data.spikesID.CellNum(s),length(data.Spikes{s}),1);
    spikeIDs(idx,3) = repmat(s,length(data.Spikes{s}),1);
end
end