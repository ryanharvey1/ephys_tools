% for laura's data
% df = readtable('F:\ClarkP30_Recordings\Analysis\hd_cell_list.csv');
% data_path = 'F:\ClarkP30_Recordings\ProcessedData\';
% save_path = 'F:\ClarkP30_Recordings\Analysis\Cell_Classification\';

df = readtable('F:\Projects\PAE_PlaceCell\analysis\swr_data\post_processed\swr_df.csv');
data_path = 'F:\Projects\PAE_PlaceCell\ProcessedData\';
save_path = 'F:\Projects\PAE_PlaceCell\analysis\cell_recruitment\';
mkdir(save_path);

% 
% df.session = df.SessionID;
df.session = df.session;

% run through each session
WaitMessage = parfor_wait(length(unique(df.session)),'Waitbar',false,'ReportInterval',1);
sessions = unique(df.session);
parfor i = 1:length(sessions)
    if exist([save_path,sessions{i},'.mat'],'file')
        continue
    end
    
%     if exist([save_path,sessions{i}],'file') % for laura's data
%         continue
%     end
    %     data = load([data_path,sessions{i},'.mat'],'avgwave','Spikes','spikesID');
    data = load([data_path,sessions{i}],'avgwave','Spikes','spikesID','frames','session_path');
    
    dat = main(data);
    
    save_data(save_path,sessions{i},dat)
    
    WaitMessage.Send;
end
WaitMessage.Destroy
%%
df = load_all(save_path);
df.waves_zscore = zscore(df.waves,0,2);


% Interneurons are selected by 2 separate criteria:
% Narrow interneuron assigned if troughToPeak <= 0.425 ms
% Wide interneuron assigned if troughToPeak > 0.425 ms and acg_tau_rise > 6 ms
% The remaining cells are assigned as Pyramidal cells.

df.cell_type(:) = {'pyr'};
df.cell_type(df.troughtoPeak <= 0.425) = {'narrow_int'};
df.cell_type(df.troughtoPeak > 0.425 & df.acg_tau_rise > 6) = {'wide_int'};


df_to_save = df;
df_to_save.waves = [];
df_to_save.waves_zscore = [];
df_to_save.acg_narrow = [];

mkdir([save_path,'processed\']);
writetable(df_to_save,[save_path,'processed\pyr_int_df.csv'])


% %
% X = [df.troughtoPeak,log10(df.spkW)];
% eva = evalclusters(X,'kmeans','CalinskiHarabasz','KList',[1:20]);
% n_clusters = eva.OptimalK;
% n_clusters = 7
% [idx,C] = kmeans(X,n_clusters,'Replicates',10,'Distance','sqeuclidean');
% 
% figure;
% h = scatterhist(df.troughtoPeak,log10(df.spkW),'Group',idx,'Kernel','on');
% grid on
% 
% 
% figure;
% for i = unique(idx)'
%     subplot(3,2,i)
%     plot(df.waves_zscore(idx==i,:)','Color',[.7,.7,.7,0.1])
%     box off
%     axis off
%     title(['cluster: ',num2str(i)])
%     darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% end
% 
% figure;
% h = scatterhist(df.troughtoPeak,log10(df.spkW),'Group',idx,'Kernel','on');
% grid on
% 
% %
% df.cell_type(:) = {'unidentified'};
% df.cell_type(df.troughtoPeak > .4 & log10(df.spkW) < 0.4) = {'pyr'};
% 
% df.cell_type(df.troughtoPeak <= .4 & log10(df.spkW) < 0.4) = {'int'};
% 
% df.cell_type(log10(df.spkW) < -0.3 & df.troughtoPeak > .4) = {'int'};
% 
% figure;
% h = scatterhist(df.troughtoPeak,log10(df.spkW),'Group',df.cell_type,'Kernel','on');
% grid on
% 
% figure;
% h = scatterhist(df.troughtoPeak,log10(df.burstIndex_Doublets),'Group',df.cell_type,'Kernel','on');
% grid on
% 
% cell_types = {'pyr','int','unidentified'};
% figure;
% for i = 1:length(cell_types)
%     subplot(3,1,i)
%     plot(df.waves_zscore(contains(df.cell_type,cell_types{i}),:)','Color',[.7,.7,.7,0.1])
%     box off
%     axis off
%     title(cell_types{i})
%     darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% end
% figure;plot(df.waves(1668,:))
% 
% %
% 
% 
% figure;
% subplot(3,1,1)
% plot(df.waves_zscore(idx==1,:)')
% subplot(3,1,2)
% plot(df.waves_zscore(idx==2,:)')
% darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% 
% 
% figure;
% h = scatterhist(df.troughtoPeak,log10(df.spkW),'Group',df.cell_type,'Kernel','on');
% grid on
% 
% figure;
% scatter(df.troughtoPeak, log10(df.acg_tau_rise),'Filled','w')
% darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% alpha .1
% 
% figure;
% scatter(df.troughtoPeak(df.cell_type=="pyr"), log10(df.acg_tau_rise(df.cell_type=="pyr")),'Filled','w')
% hold on
% scatter(df.troughtoPeak(df.cell_type=="narrow_int"),log10(df.acg_tau_rise(df.cell_type=="narrow_int")),'Filled','r')
% 
% scatter(df.troughtoPeak(df.cell_type=="wide_int"),log10(df.acg_tau_rise(df.cell_type=="wide_int")),'Filled','g')
% alpha .1
% darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% 
% 
% randn(length(df.troughToPeak(df.cell_type=="pyr")),1) * .001
% 
% figure;
% plot(df.waves_zscore(df.cell_type=="pyr",:)','Color',[.7,.7,.7,0.1])
% box off
% axis off
% darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% figure;
% plot(df.waves_zscore(df.cell_type=="narrow_int",:)','Color',[.7,.7,.7,0.1])
% box off
% axis off
% darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% figure;
% plot(df.waves_zscore(df.cell_type=="wide_int",:)','Color',[.7,.7,.7,0.1])
% box off
% axis off
% darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% %
% 
% figure;
% scatter(df.trough_to_peak(df.trough_to_peak>0.425),df.spk_wavelet(df.trough_to_peak>0.425),'filled','MarkerFaceAlpha',.2);
% hold on
% scatter(df.trough_to_peak(df.trough_to_peak<=0.425),df.spk_wavelet(df.trough_to_peak<=0.425),'filled','MarkerFaceAlpha',.2);
% 
% 
% figure;
% h = scatterhist(df.troughToPeak,df.spkW,'Group',df.cell_type,'Kernel','on');
% grid on
% 
% 
% figure;plot(df.waves(4622,:))
% figure;plot(df.waves(6141,:))
% figure;plot(df.waves(6228,:))
% 
% hold on
% plot(df.waves(5045,:) - mean(df.waves(5045,:)))
% df(5045,:)
% 
% 
% 
% 
% rough split for now
% idx = df.trough_to_peak < 0.23 ;
% df.cell_type(idx) = {'int'};
% 
% idx = df.trough_to_peak > 0.23;
% df.cell_type(idx) = {'pyr'};
% 
% df_to_save = df;
% df_to_save.waves = [];
% df_to_save.waves_zscore = [];
% mkdir([save_path,'processed\']);
% writetable(df_to_save,[save_path,'processed\pyr_int_df.csv'])
% 
% fig = figure; 
% fig.Color = [1 1 1];
% plot(df.trough_to_peak(df.trough_to_peak < 0.25& df.avg_fr > 5,1),...
%     log(df.avg_fr(df.trough_to_peak < 0.25& df.avg_fr > 5,1)),'.k')
% xlim([0 1]) 
% hold on
% plot(df.trough_to_peak(df.trough_to_peak < 0.25& df.avg_fr < 5,1),...
%     log(df.avg_fr(df.trough_to_peak < 0.25& df.avg_fr < 5,1)),'.r')
% plot(df.trough_to_peak(df.trough_to_peak < 0.25& df.avg_fr < 5,1),...
%     log(df.avg_fr(df.trough_to_peak < 0.25& df.avg_fr < 5,1)),'.r')
% plot([0.25,0.25],ylim)
% plot(xlim,[log(5),log(5)])
% WIP scratch code below
% 
% % Figure trough_to_peak vs local variation scatter
% 
% X = [df.trough_to_peak, df.local_variation,df.coeff_variation];
% X(df.short_isi > .02,:) = [];
% [coeff,score,latent,tsquared,explained,mu] = pca(X);
% figure
% scatter3(score(:,1),score(:,2),score(:,3))
% figure;
% histogram(score(:,1),100)
% 
% 
% % Z-scored waveforms 
% figure;
% idx = df.avg_fr > 5 & df.trough_to_peak < 0.23;
% plot(df.waves_zscore(idx,:)','Color',[.7,.7,.7,0.1])
% box off
% axis off
% darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% 
% figure;
% idx = df.avg_fr < 5 & df.trough_to_peak > 0.23;
% plot(df.waves_zscore(idx,:)','Color',[.7,.7,.7,0.1])
% box off
% axis off
% darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% 
% figure;
% idx = df.avg_fr < 5 & df.trough_to_peak < 0.23;
% plot(df.waves_zscore(idx,:)','Color',[.7,.7,.7,0.1])
% box off
% axis off
% %
% X = [df.trough_to_peak, df.spk_wavelet,log(df.coeff_variation)];
% figure; plotmatrix(X)
% 
% [coeff,score,latent,tsquared,explained,mu] = pca(X);
% figure
% scatter(score(:,1),score(:,2))
% figure;
% histogram(score(:,1),100)
% 
% eva = evalclusters(score(:,1:2),'kmeans','CalinskiHarabasz','KList',[1:6]);
% eva.OptimalK
% [idx,C] = kmeans(score(:,1:2),4,'Replicates',10,'Distance','sqeuclidean');
% 
% 
% figure;
% subplot(3,1,1)
% plot(df.waves_zscore(idx==1,:)')
% subplot(3,1,2)
% plot(df.waves_zscore(idx==2,:)')
% darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% 
% figure
% scatter(score(idx==1,1),score(idx==1,2),5,'filled')
% hold on
% scatter(score(idx==2,1),score(idx==2,2),5,'filled')
% scatter(score(idx==3,1),score(idx==3,2),5,'filled')
% 
% darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% 
% figure;
% scatter(df.local_variation(idx==1),df.trough_to_peak(idx==1),5,'filled')
% hold on
% scatter(df.local_variation(idx==2),df.trough_to_peak(idx==2),5,'filled')
% 
% %
% 
% df.waves_zscore = zscore(df.waves,0,2);
% idx = df.trough_to_peak > 0.25;
% figure;
% plot(df.waves_zscore(idx,:)')
% 
% idx = df.trough_to_peak < 0.25;
% figure;
% plot(df.waves_zscore(idx,:)')
% %
% 
% X = [df.peak_to_trough, df.trough_to_peak, df.spk_wavelet, df.local_variation,log(df.coeff_variation)];
% X = normalize(X,'norm');
% 
% [coeff,score,latent,tsquared,explained,mu] = pca(X);
% rng(3)
% k = 2;
% d = 500;
% x1 = linspace(min(score(:,1)) - 2,max(score(:,1)) + 2,d);
% x2 = linspace(min(score(:,2)) - 2,max(score(:,2)) + 2,d);
% [x1grid,x2grid] = meshgrid(x1,x2);
% X0 = [x1grid(:) x2grid(:)];
% threshold = sqrt(chi2inv(0.99,2));
% options = statset('MaxIter',1000); % Increase number of EM iterations
% AIC = zeros(2,1);
% gmfit = cell(2,1);
% 
% figure;
% for k=1:2
%     gmfit{k} = fitgmdist(score(:,1:2),k,'CovarianceType','diagonal',...
%         'SharedCovariance',true,'Options',options);
%     clusterX = cluster(gmfit{k},score(:,1:2));
%     mahalDist = mahal(gmfit{k},X0);
%     subplot(1,2,k);
%     h1 = gscatter(score(:,1),score(:,2),clusterX);
%     hold on;
%     for m = 1:k
%         idx = mahalDist(:,m)<=threshold;
%         Color = h1(m).Color*0.75 + -0.5*(h1(m).Color - 1);
%         h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
%         uistack(h2,'bottom');
%     end
%     plot(gmfit{k}.mu(:,1),gmfit{k}.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
%     title(sprintf('Sigma is %s, SharedCovariance = %s',...
%         'diagonal',true),'FontSize',8)
%     legend(h1,{'1','2','3'});
%     hold off
%     AIC(k)= GMModels{k}.AIC;
%     AIC(k) = gmfit{k}.AIC ;
% end
% darkBackground(gcf,[0.2 0.2 0.2],[1 1 1])
% [minAIC,numComponents] = min(AIC);
% numComponents
% FisrtModel = GMModels{1}
% 
% 
% %
% X = [df.peak_to_trough, df.trough_to_peak, df.spk_wavelet, df.local_variation,log(df.avg_fr)];
% X = normalize(X,'norm');
% 
% [reduction, umap, clusterIdentifiers] = run_umap(X);
% [reduction, umap, clusterIdentifiers] = run_umap(df.waves_zscore);
% 
% clusters = unique(clusterIdentifiers);
% p = ceil(sqrt(length(clusters)));
% figure;
% for i = double(clusters)
%     subplot(p,p,i+1)
%     plot(df.waves_zscore(clusterIdentifiers == i,:)','Color',[.7,.7,.7,0.1])
%     box off
%     axis off
% end
% darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
% 
% idx = clusterIdentifiers==34 | clusterIdentifiers==35;
% figure;
% scatter(df.trough_to_peak,df.avg_fr,5,'k','filled')
% hold on
% scatter(df.trough_to_peak(idx),df.avg_fr(idx),5,'r','filled')
% 
% df.session(clusterIdentifiers==5)
% %
% meas = X;
% Y = tsne(X,'Algorithm','barneshut','Distance','mahalanobis');
% subplot(2,2,1)
% figure
% scatter(Y(:,1),Y(:,2))
% title('Mahalanobis')
% 
% rng('default') % for fair comparison
% Y = tsne(X,'Algorithm','barneshut','Distance','cosine');
% subplot(2,2,2)
% figure
% scatter(Y(:,1),Y(:,2))
% title('Cosine')
% 
% rng('default') % for fair comparison
% Y = tsne(X,'Algorithm','barneshut','Distance','chebychev');
% subplot(2,2,3)
% scatter(Y(:,1),Y(:,2))
% title('Chebychev')
% 
% rng('default') % for fair comparison
% Y = tsne(X,'Algorithm','barneshut','Distance','euclidean','Standardize',true);
% subplot(2,2,4)
% scatter(Y(:,1),Y(:,2))
% title('Euclidean')
%%
% local functions below

function dat=main(data)

session_info = LoadParameters(data.session_path);
data.fs = session_info.rates.wideband;

% put the waveforms into a matrix and make sure they are 32 samples long (1ms)
waves = get_waveform(data);
% flip waveforms so the trough is pointing down
% (neuralynx flips signal by default. Also, waveforms can be flipped across
% a dipole)

waveform_metrics = calc_waveform_metrics(data);

% [waves,flipped] = get_flipped_waveform(waveform_metrics.waves);
% get waveform peak to trough and trough to peak times
% [peakToTrough,troughToPeak] = get_waveform_delay_times(waveform_metrics.waves);
% get spike duration via wavelet transform
spkW = get_spike_wavelet(waves,waveform_metrics.troughtoPeak,data.fs);
% get short isi,local variation (DOI: 10.1162/089976603322518759),
% coefficient of variation.
[short_isi,lv,cv] = isi_metrics(data.Spikes);
% get average firing rate
avg_fr = get_average_fr(data.Spikes);
% get waveform asymmetry
asymmetry = get_waveform_asymmetry(waveform_metrics.waves');
% get number of spikes
n_spikes = cellfun(@numel, data.Spikes)';
% get acg metricts 
acg_metrics = calc_ACG_metrics(data);
dat = fit_ACG(acg_metrics.acg_narrow);

for names = fields(acg_metrics)'
   dat.(names{1}) =  acg_metrics.(names{1});
end

for names = fields(waveform_metrics)'
   dat.(names{1}) =  waveform_metrics.(names{1});
end

% dat.waves = waveform_metrics.waves';
% dat.flipped = flipped';
% dat.peakToTrough = peakToTrough';
% dat.troughToPeak = troughToPeak';
dat.spkW = spkW';
dat.short_isi = short_isi; 
dat.lv = lv;
dat.cv = cv;
dat.id = data.spikesID;
dat.avg_fr = avg_fr;
dat.asymmetry = asymmetry;
dat.n_spikes = n_spikes;
end


function waveform_metrics = calc_waveform_metrics(data)
% Extracts waveform metrics
% 
% INPUT:
% waveforms structure
% 
% OUTPUT:
% waveform_metrics
% Metrics for the waveforms: peaktoTrough, troughtoPeak,
% derivative_TroughtoPeak, peakA, peakB, ab_ratio, trough
  
% By Peter Petersen
% petersen.peter@gmail.com

% format data
for i = 1:length(data.avgwave)
    [~,max_idx] = max(max(abs(data.avgwave{i}),[],2));
    wave = data.avgwave{i}(max_idx,:);
    waveforms.filtWaveform{i} = interp1(linspace(1,32,length(wave)),wave,1:32);
    [~,loc(i)] = max(abs(waveforms.filtWaveform{i}));
end
sr_in = data.fs;
waveforms.timeWaveform{1} = ((1:32) - mode(loc)) / sr_in * 1000;


filtWaveform = waveforms.filtWaveform;
timeWaveform = waveforms.timeWaveform{1};
timeWaveform_span = length(timeWaveform) * mean(diff(timeWaveform));
% sr = 1/mean(diff(timeWaveform))*1000;
oversampling = ceil(100000/sr_in);
sr = oversampling * sr_in;
timeWaveform = interp1(timeWaveform,timeWaveform,timeWaveform(1):mean(diff(timeWaveform))/oversampling:timeWaveform(end),'spline');
zero_idx = find(timeWaveform>=0,1);
trough_interval = [find(timeWaveform>=-0.25,1),find(timeWaveform>=0.25,1)]; % -10/40:10/40 => 
trough_interval_1stDerivative = [find(timeWaveform>=-0.50,1),find(timeWaveform>=0.25,1)]; % -20/40:10/40 => 
trough_interval_2stDerivative = [find(timeWaveform>=-0.625,1),find(timeWaveform>=-0.125,1)]; % 7:27 => -25/40:-5/40 => -0.625:-0.125
idx_45 = find(timeWaveform>=0.325,1);
idx_49 = find(timeWaveform>=0.425,1);
idx_3 = find(timeWaveform>= -0.725,1);
idx_6 = find(timeWaveform>= -0.650,1);
idx_9 = find(timeWaveform>= -0.575,1);
idx_end4 = find(timeWaveform>= 0.700,1);
idx_end2 = find(timeWaveform>= 0.750,1);
waveform_metrics = [];
m = 0;
n = 0;
wave = [];
wave_diff = [];
wave_diff2 = [];
wave_cut = [];
wave_align = [];
wave_diff_cut = [];

t_before = [];
t_after = [];
peakA = [];
peakB = [];
trough = [];
width1 = [];
t_after_diff = [];
Itest = [];
Itest2 = [];

for m = 1:length(filtWaveform)
    if ~any(isnan(filtWaveform{m}))
        wave = interp1(waveforms.timeWaveform{1},zscore(filtWaveform{m}),timeWaveform(1):mean(diff(timeWaveform)):timeWaveform(end),'spline');
        waveform_metrics.polarity(m) = mean(wave(trough_interval(1):trough_interval(2))) - mean(wave([1:trough_interval(1),trough_interval(2):end]));
        if waveform_metrics.polarity(m) > 0
            wave = -wave;
        end
        wave_diff{m} = diff(wave);
        wave_diff2{m} = diff(wave,2);
        [MIN2,I2] = min(wave(trough_interval(1):trough_interval(2))); % trough_interval
        [MIN2_diff,I2_diff] = min(wave_diff{m}(trough_interval_1stDerivative(1):trough_interval_1stDerivative(2))); % 1st deriv
        [MIN2_diff2,I2_diff2] = min(wave_diff2{m}(trough_interval_2stDerivative(1):trough_interval_2stDerivative(2))); % 2nd deriv, trough_interval_2stDerivative
        [MAX3,I3] = max(wave(1:I2+trough_interval(1)-1));
        [MAX4,I4] = max(wave(I2+trough_interval(1):end));
        [MAX4_diff,I4_diff] = max(wave_diff{m}(I2_diff+trough_interval_1stDerivative(1):end));
        t_before(m) = I2+trough_interval(1)-1-I3;
        t_after(m) = I4;
        t_after_diff(m) = I4_diff;
        peakA(m) = MAX3;
        peakB(m) = MAX4;
        trough(m) = MIN2;
        Itest(m) = I2;
        Itest2(m) = length(wave)-(idx_45+I2);
        indexes = I3:I2+I4+trough_interval(1)-1;
        
        waves(:,m) = wave;
    else
        t_before(m) = nan;
        t_after(m) = nan;
        t_after_diff(m) = nan;
        peakA(m) = nan;
        peakB(m) = nan;
        trough(m) = nan;
    end
end
waveform_metrics.peaktoTrough = t_before/sr*1000;
waveform_metrics.troughtoPeak = t_after/sr*1000;
waveform_metrics.derivative_TroughtoPeak = t_after_diff/sr*1000;
waveform_metrics.peakA = peakA;
waveform_metrics.peakB = peakB;
waveform_metrics.ab_ratio = (peakB-peakA)./(peakA+peakB);
waveform_metrics.trough = trough;
waveform_metrics.waves = waves;
end

% function [acg_metrics] = calc_acg_metrics(temp)
% session_info = LoadParameters(temp.session_path);
% 
% [ccg,time] = CCG(temp.Spikes,[],'binSize',0.001,'duration',1,'norm','rate','Fs',session_info.rates.wideband);
% [ccg2,time2] = CCG(temp.Spikes,[],'binSize',0.0005,'duration',0.100,'norm','rate','Fs',session_info.rates.wideband);
% max_lags = (length(time)-1)/2;
% max_lags2 = (length(time2)-1)/2;
% for ii = 1:size(ccg,2)
%     xc = ccg(:,ii,ii);
%     xc2 = ccg2(:,ii,ii);
%     ThetaModulationIndex(ii) = (mean(xc(max_lags+1+50:max_lags+1+70))-...
%         mean(xc(max_lags+1+100:max_lags+1+140)))/...
%         (mean(xc(max_lags+1+50:max_lags+1+70))+...
%         mean(xc(max_lags+1+100:max_lags+1+140)));
%     BurstIndex_Royer2012(ii,1) = sum(xc(max_lags+1+3:max_lags+1+5))/...
%         mean(xc(max_lags+1+200:max_lags+1+300));
%     BurstIndex_Doublets(ii,1) = max(xc2(max_lags2+1+5:max_lags2+1+16))/...
%         mean(xc2(max_lags2+1+16:max_lags2+1+23));
%     ACG(:,ii) = xc;
%     ACG2(:,ii) = ccg2(:,ii,ii);
% end
% acg_metrics.acg = ACG;
% acg_metrics.acg2 = ACG2;
% acg_metrics.thetaModulationIndex = ThetaModulationIndex;
% acg_metrics.burstIndex_Royer2012 = BurstIndex_Royer2012;
% acg_metrics.burstIndex_Doublets = BurstIndex_Doublets;
% acg_metrics.ccg = ccg2;
% acg_metrics.ccg_time = time2;
% end

function acg_metrics = calc_ACG_metrics(data)%(clustering_path,sr,TimeRestriction)
% Two autocorrelograms are calculated:  narrow (100ms, 0.5ms bins) and wide (1s, 1ms bins) using the CCG function (for speed)
%
% Further three metrics are derived from these:
%
% Theta modulation index:
%    Computed as the difference between the theta modulation trough (defined as mean of autocorrelogram bins 50-70 ms)
%    and the theta modulation peak (mean of autocorrelogram  bins 100-140ms) over their sum
%    Originally defined in Cacucci et al., JNeuro 2004
%
% BurstIndex_Doublets:
%    max bin count from 2.5-8ms normalized by the average number of spikes in the 8-11.5ms bins
%
% BurstIndex_Royer2012:
%    Burst index is determined by calculating the average number of spikes in the 3-5 ms bins of the spike
%    autocorrelogram divided by the average number of spikes in the 200-300 ms bins.
%    Metrics introduced in Royer et al. Nature Neuroscience 2012, and adjusted in Senzai & Buzsaki, Neuron 2017.

% By Peter Petersen
% petersen.peter@gmail.com
% Last edited: 06-10-2020

spikes.numcells = length(data.Spikes);
spikes.times = data.Spikes;
sr = data.fs;

ThetaModulationIndex = nan(1,spikes.numcells);
BurstIndex_Royer2012 = nan(1,spikes.numcells);
BurstIndex_Doublets = nan(1,spikes.numcells);


cell_indexes = 1:spikes.numcells;

bins_wide = 500;
acg_wide = zeros(bins_wide*2+1,numel(spikes.times));
bins_narrow = 100;
acg_narrow = zeros(bins_narrow*2+1,numel(spikes.times));
% disp('Calculating narrow ACGs (100ms, 0.5ms bins) and wide ACGs (1s, 1ms bins)')
% tic
for i = cell_indexes
    acg_wide(:,i) = CCG(spikes.times{i},ones(size(spikes.times{i})),'binSize',0.001,'duration',1,'norm','rate','Fs',1/sr);
    acg_narrow(:,i) = CCG(spikes.times{i},ones(size(spikes.times{i})),'binSize',0.0005,'duration',0.100,'norm','rate','Fs',1/sr);
    % Metrics from narrow
    BurstIndex_Doublets(i) = max(acg_narrow(bins_narrow+1+5:bins_narrow+1+16,i))/mean(acg_narrow(bins_narrow+1+16:bins_narrow+1+23,i));
    % Metrics from wide
    ThetaModulationIndex(i) = (mean(acg_wide(bins_wide+1+100:bins_wide+1+140,i)) - mean(acg_wide(bins_wide+1+50:bins_wide+1+70,i)))/(mean(acg_wide(bins_wide+1+50:bins_wide+1+70,i))+mean(acg_wide(bins_wide+1+100:bins_wide+1+140,i)));
    BurstIndex_Royer2012(i) = mean(acg_wide(bins_wide+1+3:bins_wide+1+5,i))/mean(acg_wide(bins_wide+1+200:bins_wide+1+300,i));
end
% toc

% figure, subplot(3,1,1)
% histogram(ThetaModulationIndex,40),xlabel('Theta modulation index'), ylabel('Count')
% subplot(3,1,2)
% histogram(BurstIndex_Royer2012,40),xlabel('BurstIndex Royer2012'), ylabel('Count')
% subplot(3,1,3)
% histogram(BurstIndex_Doublets,40),xlabel('BurstIndex Doublets'), ylabel('Count')

acg_metrics.acg_wide = acg_wide;
acg_metrics.acg_narrow = acg_narrow;
acg_metrics.thetaModulationIndex = ThetaModulationIndex;
acg_metrics.burstIndex_Royer2012 = BurstIndex_Royer2012;
acg_metrics.burstIndex_Doublets = BurstIndex_Doublets;
end

function fit_params_out = fit_ACG(acg_narrow)
% This function is part of CellExplorer
% Fits a tripple exponential to the autocorrelogram with 0.5ms bins from -50ms -> 50ms 

% By Peter Petersen
% Last edited: 24-08-22020;

acg_narrow(100:102) = 0; % Sets the time-zero bin to zero (-0.5ms -> 0.5ms)
fit_params = nan(8,size(acg_narrow,2));
rsquare = nan(1,size(acg_narrow,2));
plotf0 = zeros(100,size(acg_narrow,2));
offset = 101;
x = ([1:100]/2)';
g = fittype('max(c*(exp(-(x-f)/a)-d*exp(-(x-f)/b))+h*exp(-(x-f)/g)+e,0)',...
    'dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c','d','e','f','g','h'});

% Turning off interation limit warning
warning('off','stats:nlinfit:IterationLimitExceeded')

% Fitting ACGs in parfor
gcp;
parfor j = 1:size(acg_narrow,2)
    if ~any(isnan(acg_narrow(:,j)))
        [f0,gof] = fit(x,acg_narrow(x*2+offset,j),g,...
            'StartPoint',[20, 1, 30, 2, 0, 5, 1.5,2],...
            'Lower',[1, 0.1, 0, 0, -30, 0,0.1,0],...
            'Upper',[500, 50, 500, 15, 50, 20,5,100]);
        plotf0(:,j) = f0(x);
        fit_params(:,j) = coeffvalues(f0);
        rsquare(j) = gof.rsquare;
    end
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

% function fit_params_out = fit_ACG(acg2)
% figs = false;
% 
% fit_params = [];
% rsquare = [];
% offset = 101;
% 
% g = fittype('max(c*(exp(-(x-f)/a)-d*exp(-(x-f)/b))+h*exp(-(x-f)/g)+e,0)',...
%     'dependent',{'y'},'independent',{'x'},'coefficients',{'a','b','c','d','e','f','g','h'});
% 
% acg_mat =[]; paut=[];
% jj = 1;
% for j = 1:size(acg2,2)
%     x = ([1:100]/2)';
%     y = acg2(x*2+offset,j); % /max(acg2(x+offset,j))
%     
%     [f0,gof] = fit(x,y,g,'StartPoint',[20, 1, 30, 2, 0, 5, 1.5,2],...
%         'Lower',[1, 0.1, 0, 0, -30, 0,0.1,0],...
%         'Upper',[500, 50, 500, 15, 50, 20,5,100]);
%     %     f0 = fit(x,y,'exp2','StartPoint',[1, -0.015, -1, -1]);
%     fit_params(:,j) = coeffvalues(f0);
%     %     xx = linspace(0,48,100);
%     rsquare(j) = gof.rsquare;
% end
% fit_params_out.acg_tau_decay = fit_params(1,:);
% fit_params_out.acg_tau_rise = fit_params(2,:);
% fit_params_out.acg_c = fit_params(3,:);
% fit_params_out.acg_d = fit_params(4,:);
% fit_params_out.acg_asymptote = fit_params(5,:);
% fit_params_out.acg_refrac = fit_params(6,:);
% fit_params_out.acg_tau_burst = fit_params(7,:);
% fit_params_out.acg_h = fit_params(8,:);
% fit_params_out.acg_fit_rsquare = rsquare;
% end

function waves = get_waveform(temp)
for w = 1:length(temp.avgwave)
    temp.avgwave{w} = interp1(linspace(1,32,length(temp.avgwave{w})),...
        vertcat(temp.avgwave{w})',1:32)';
    [~,I] = max(max(abs(temp.avgwave{w}),[],2));
    waves(w,:) = temp.avgwave{w}(I,:);
end
end

function [waveforms,flipped] = get_flipped_waveform(waveforms)
for a = 1:size(waveforms,1)
    thiswave = waveforms(a,:);
    [~,minpos] = min(thiswave); %find the trough
    minpos = minpos(1);
    
    [~,maxpos] = max(thiswave); %find the largest peak
    maxpos = maxpos(1);
    
    if abs(thiswave(maxpos)) > abs(thiswave(minpos))
        %flip waveform so peak becomes trough - for waveforms that have large trough preceeding peak.
        thiswave = thiswave*-1;
        [~,minpos] = min(thiswave); %recalculate trough
        minpos = minpos(1);
        [~,maxpos] = max(thiswave); %recalculate peak
        maxpos = maxpos(1);
        
        % save flipped waveform
        waveforms(a,:) = thiswave;
        
        % save flipped indicator
        flipped(a,1) = true;
    else
        % save flipped indicator
        flipped(a,1) = false;
    end
end
end

function [peakToTrough,troughToPeak]=get_waveform_delay_times(waveforms)

for a = 1:size(waveforms,1)
    
    thiswave = waveforms(a,:);
    
    %find the trough
    [~,minpos] = min(thiswave);
    minpos = minpos(1);
    
    %find the largest peak
    [~,maxpos] = max(thiswave);
    maxpos = maxpos(1);
    
    % Find time from min to left %finds peak before trough
    [max_left,~] = max(thiswave(1,1:minpos-1));
    if isempty(max_left)
        maxpos_left = minpos;
    else
        maxpos_left = find(thiswave == max_left(1));
    end
    
    peakToTrough(a,:) = (minpos - maxpos_left(1))/length(thiswave); % in ms
    
    % Find time from min to right
    [max_right,~] = max(thiswave(1,minpos+1:end)); %finds peak after trough
    if isempty(max_right)
        max_right = minpos;
    else
        maxpos_right = find(thiswave == max_right(1));
    end
    
    troughToPeak(a,:) = (maxpos_right(1) -  minpos)/length(thiswave); % in ms
    
    
    %         figure;
    %         plot(thiswave);
    %         hold on;
    %         plot(minpos,thiswave(minpos),'*b'); hold on;
    %         plot(maxpos,thiswave(maxpos),'*r'); hold on;
    %         plot(maxpos_right,thiswave(maxpos_right),'*k');
    %         plot(maxpos_left,thiswave(maxpos_left),'*c');
    %         title(['peak2trough =  ', num2str(peakToTrough(a,:)),'  trough2peak =  ', num2str(troughToPeak(a,:))])
end
end


function spkW = get_spike_wavelet(waveforms,peakToTrough,fs)

for a = 1:size(waveforms,1)
    w = waveforms(a,:)';
    if isnan(peakToTrough(a))
        spkW(a,1) = NaN;
        continue
    end
    
    w = [w(1)*ones(1000,1);w;w(end)*ones(1000,1)];
    [wave, f, t] = getWavelet(w,fs,300,8000,128);
    %We consider only the central portion of the wavelet because we
    %haven't filtered it before hand (e.g. with a Hanning window)
    wave = wave(:,int16(length(t)/4):3*int16(length(t)/4));
    %Where is the max frequency?
    [maxPow, ix] = max(wave);
    [~, mix] = max(maxPow);
    ix = ix(mix);
    spkW(a,1) = 1000/f(ix);
end
end

function [short_isi,lv,cv] = isi_metrics(spikes)
for i = 1:size(spikes,1)
    ts = spikes{i};
    T = diff(ts);
    lag = .002;
    short_isi(i) = sum(T < lag)/length(T);
    T = T(T > lag);
    lv(i) = (1/(length(T)-1)) * nansum((3*(diff(T)).^2) ./ ((T(1:end-1) + T(2:end)).^2));
    cv(i) = std(T)/mean(T);
end
end

function avg_fr = get_average_fr(spikes)
all_spikes = vertcat(spikes{:});
duration = max(all_spikes) - min(all_spikes);
for i = 1:size(spikes,1)
    avg_fr(i) = length(spikes{i}) / duration;
end
end

function asymmetry = get_waveform_asymmetry(waveforms)
for a = 1:size(waveforms,1)
    
    thiswave = zscore(waveforms(a,:));
    
%     try
%         [pks,locs,widths,proms] = findpeaks(thiswave,...
%             'MinPeakHeight',max(thiswave)*.6,'MinPeakDistance',20);
%         [~,idx] = sort(locs,'descend');
%         pks = pks(idx);
%         locs = locs(idx);
%         widths = widths(idx);
%         proms = proms(idx);
%         
%         peak_area = abs(proms(1) * widths(1));
%         peak_val = abs(pks(1));
%         
%         [pks,locs,widths,proms] = findpeaks(thiswave*-1,...
%             'MinPeakHeight',max(thiswave*-1)*.6,'MinPeakDistance',20);
%         [~,idx] = sort(locs,'descend');
%         pks = pks(idx);
%         locs = locs(idx);
%         widths = widths(idx);
%         proms = proms(idx);
%         
%         trough_area = abs(proms(1) * widths(1));
%         trough_val = abs(pks(1));
%         
%         asymmetry(a) = (trough_val / peak_val) * (trough_area/peak_area);
%     catch
        
        [trough_val,trough_pos] = min(thiswave);
        
        [peak_val,~] = max(thiswave(trough_pos:end));
        
        asymmetry(a) = (abs(trough_val) / abs(peak_val)) *...
            (trapz(abs(thiswave(thiswave < 0))) / trapz(abs(thiswave(thiswave > 0))));
%     end
    
end
end

% function save_data(save_path,session,dat)
% save([save_path,session,'.mat'],'dat')
% end


function save_data(save_path,session,dat)
save([save_path,session],'dat')
end

function df = load_all(save_path)
df = table();
file = dir([save_path,'*.mat']);
for i = 1:length(file)
    load([save_path,file(i).name],'dat');
    if ~isstruct(dat)
        continue
    end
    df=[df;construct_df(dat,file(i).name)];
end
end

function df=construct_df(dat,session)
df = table();
% set up id
df.session = repmat({extractBefore(session,'.mat')},length(dat.lv),1);
df.tetrode = dat.id.TetrodeNum;
df.cell = dat.id.CellNum;
% extract metrics

% order dat so matrices are at the end
[~,idx]= sort(structfun(@(x) size(x,1), dat));
dat = orderfields(dat,idx);

% pull out each field and add to df
for names = fields(dat)'
    if contains(names,'id')
        continue
    end
   df.(names{1}) =  dat.(names{1})';
end

% df.flipped = dat.flipped;
% df.peak_to_trough = dat.peakToTrough;
% df.trough_to_peak = dat.troughToPeak;
% df.spk_wavelet = dat.spkW;
% df.local_variation = dat.lv';
% df.coeff_variation = dat.cv'; 
% df.short_isi = dat.short_isi'; 
% df.avg_fr = dat.avg_fr';
% df.asymmetry = dat.asymmetry';
% df.n_spikes = dat.n_spikes';
% % add waveforms
% df.waves = dat.waves;
end
