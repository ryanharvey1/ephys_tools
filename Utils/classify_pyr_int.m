
df = readtable('F:\ClarkP30_Recordings\Analysis\hd_cell_list.csv');
data_path = 'F:\ClarkP30_Recordings\ProcessedData\';
save_path = 'F:\ClarkP30_Recordings\Analysis\Cell_Classification\';

% 
df.session = df.SessionID;
% run through each session
WaitMessage = parfor_wait(length(unique(df.session)),'Waitbar',false,'ReportInterval',1);
sessions = unique(df.session);
parfor i = 1:length(sessions)
    %     if exist([save_path,sessions{i},'.mat'],'file')
    %         continue
    %     end
    
    if exist([save_path,sessions{i}],'file')
        continue
    end
    %     data = load([data_path,sessions{i},'.mat'],'avgwave','Spikes','spikesID');
    data = load([data_path,sessions{i}],'avgwave','Spikes','spikesID','frames');
    
    dat = main(data);
    
    save_data(save_path,sessions{i},dat)
    
    WaitMessage.Send;
end
WaitMessage.Destroy
%%
df = load_all(save_path);
df.waves_zscore = zscore(df.waves,0,2);

% rough split for now
idx = df.trough_to_peak < 0.23 ;
df.cell_type(idx) = {'int'};

idx = df.trough_to_peak > 0.23;
df.cell_type(idx) = {'pyr'};

df_to_save = df;
df_to_save.waves = [];
df_to_save.waves_zscore = [];
mkdir([save_path,'processed\']);
writetable(df_to_save,[save_path,'processed\pyr_int_df.csv'])

fig = figure; 
fig.Color = [1 1 1];
plot(df.trough_to_peak(df.trough_to_peak < 0.25& df.avg_fr > 5,1),...
    log(df.avg_fr(df.trough_to_peak < 0.25& df.avg_fr > 5,1)),'.k')
xlim([0 1]) 
hold on
plot(df.trough_to_peak(df.trough_to_peak < 0.25& df.avg_fr < 5,1),...
    log(df.avg_fr(df.trough_to_peak < 0.25& df.avg_fr < 5,1)),'.r')
plot(df.trough_to_peak(df.trough_to_peak < 0.25& df.avg_fr < 5,1),...
    log(df.avg_fr(df.trough_to_peak < 0.25& df.avg_fr < 5,1)),'.r')
plot([0.25,0.25],ylim)
plot(xlim,[log(5),log(5)])
% WIP scratch code below

%% Figure trough_to_peak vs local variation scatter

X = [df.trough_to_peak, df.local_variation,df.coeff_variation];
X(df.short_isi > .02,:) = [];
[coeff,score,latent,tsquared,explained,mu] = pca(X);
figure
scatter3(score(:,1),score(:,2),score(:,3))
figure;
histogram(score(:,1),100)


%% Z-scored waveforms 
figure;
idx = df.avg_fr > 5 & df.trough_to_peak < 0.23;
plot(df.waves_zscore(idx,:)','Color',[.7,.7,.7,0.1])
box off
axis off
darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])

figure;
idx = df.avg_fr < 5 & df.trough_to_peak > 0.23;
plot(df.waves_zscore(idx,:)','Color',[.7,.7,.7,0.1])
box off
axis off
darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])

figure;
idx = df.avg_fr < 5 & df.trough_to_peak < 0.23;
plot(df.waves_zscore(idx,:)','Color',[.7,.7,.7,0.1])
box off
axis off
%%
X = [df.trough_to_peak, df.spk_wavelet,log(df.coeff_variation)];
figure; plotmatrix(X)

[coeff,score,latent,tsquared,explained,mu] = pca(X);
figure
scatter(score(:,1),score(:,2))
figure;
histogram(score(:,1),100)

eva = evalclusters(score(:,1:2),'kmeans','CalinskiHarabasz','KList',[1:6]);
% eva.OptimalK
[idx,C] = kmeans(score(:,1:2),4,'Replicates',10,'Distance','sqeuclidean');


figure;
subplot(3,1,1)
plot(df.waves_zscore(idx==1,:)')
subplot(3,1,2)
plot(df.waves_zscore(idx==2,:)')
darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])

figure
scatter(score(idx==1,1),score(idx==1,2),5,'filled')
hold on
scatter(score(idx==2,1),score(idx==2,2),5,'filled')
scatter(score(idx==3,1),score(idx==3,2),5,'filled')

darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])

figure;
scatter(df.local_variation(idx==1),df.trough_to_peak(idx==1),5,'filled')
hold on
scatter(df.local_variation(idx==2),df.trough_to_peak(idx==2),5,'filled')

%%

df.waves_zscore = zscore(df.waves,0,2);
idx = df.trough_to_peak > 0.25;
figure;
plot(df.waves_zscore(idx,:)')

idx = df.trough_to_peak < 0.25;
figure;
plot(df.waves_zscore(idx,:)')
%%

X = [df.peak_to_trough, df.trough_to_peak, df.spk_wavelet, df.local_variation,log(df.coeff_variation)];
X = normalize(X,'norm');

[coeff,score,latent,tsquared,explained,mu] = pca(X);
rng(3)
k = 2;
d = 500;
x1 = linspace(min(score(:,1)) - 2,max(score(:,1)) + 2,d);
x2 = linspace(min(score(:,2)) - 2,max(score(:,2)) + 2,d);
[x1grid,x2grid] = meshgrid(x1,x2);
X0 = [x1grid(:) x2grid(:)];
threshold = sqrt(chi2inv(0.99,2));
options = statset('MaxIter',1000); % Increase number of EM iterations
AIC = zeros(2,1);
gmfit = cell(2,1);

figure;
for k=1:2
    gmfit{k} = fitgmdist(score(:,1:2),k,'CovarianceType','diagonal',...
        'SharedCovariance',true,'Options',options);
    clusterX = cluster(gmfit{k},score(:,1:2));
    mahalDist = mahal(gmfit{k},X0);
    subplot(1,2,k);
    h1 = gscatter(score(:,1),score(:,2),clusterX);
    hold on;
    for m = 1:k
        idx = mahalDist(:,m)<=threshold;
        Color = h1(m).Color*0.75 + -0.5*(h1(m).Color - 1);
        h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
        uistack(h2,'bottom');
    end
    plot(gmfit{k}.mu(:,1),gmfit{k}.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
    title(sprintf('Sigma is %s, SharedCovariance = %s',...
        'diagonal',true),'FontSize',8)
    legend(h1,{'1','2','3'});
    hold off
    % AIC(k)= GMModels{k}.AIC;
    AIC(k) = gmfit{k}.AIC ;
end
darkBackground(gcf,[0.2 0.2 0.2],[1 1 1])
[minAIC,numComponents] = min(AIC);
numComponents
% FisrtModel = GMModels{1}


%%
X = [df.peak_to_trough, df.trough_to_peak, df.spk_wavelet, df.local_variation,log(df.avg_fr)];
X = normalize(X,'norm');

% [reduction, umap, clusterIdentifiers] = run_umap(X);
[reduction, umap, clusterIdentifiers] = run_umap(df.waves_zscore);

clusters = unique(clusterIdentifiers);
p = ceil(sqrt(length(clusters)));
figure;
for i = double(clusters)
    subplot(p,p,i+1)
    plot(df.waves_zscore(clusterIdentifiers == i,:)','Color',[.7,.7,.7,0.1])
    box off
    axis off
end
darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])

idx = clusterIdentifiers==34 | clusterIdentifiers==35;
figure;
scatter(df.trough_to_peak,df.avg_fr,5,'k','filled')
hold on
scatter(df.trough_to_peak(idx),df.avg_fr(idx),5,'r','filled')

df.session(clusterIdentifiers==5)
%%
% meas = X;
Y = tsne(X,'Algorithm','barneshut','Distance','mahalanobis');
subplot(2,2,1)
figure
scatter(Y(:,1),Y(:,2))
title('Mahalanobis')

rng('default') % for fair comparison
Y = tsne(X,'Algorithm','barneshut','Distance','cosine');
subplot(2,2,2)
figure
scatter(Y(:,1),Y(:,2))
title('Cosine')

rng('default') % for fair comparison
Y = tsne(X,'Algorithm','barneshut','Distance','chebychev');
subplot(2,2,3)
scatter(Y(:,1),Y(:,2))
title('Chebychev')

rng('default') % for fair comparison
Y = tsne(X,'Algorithm','barneshut','Distance','euclidean','Standardize',true);
subplot(2,2,4)
scatter(Y(:,1),Y(:,2))
title('Euclidean')
%%
% local functions below

function dat=main(data)
% put the waveforms into a matrix and make sure they are 64 samples long (2ms)
waves = get_waveform(data);
% flip waveforms so the trough is pointing down
% (neuralynx flips signal by default. Also, waveforms can be flipped across
% a dipole)
[waves,flipped] = get_flipped_waveform(waves);
% get waveform peak to trough and trough to peak times
[peakToTrough,troughToPeak] = get_waveform_delay_times(waves);
% get spike duration via wavelet transform
spkW = get_spike_wavelet(waves,peakToTrough);
% get short isi,local variation (DOI: 10.1162/089976603322518759),
% coefficient of variation.
 [short_isi,lv,cv] = isi_metrics(data.Spikes);
% get average firing rate
avg_fr = get_average_fr(data.Spikes);
% get waveform asymmetry
asymmetry = get_waveform_asymmetry(waves);



dat.waves = waves;
dat.flipped = flipped;
dat.peakToTrough = peakToTrough;
dat.troughToPeak = troughToPeak;
dat.spkW = spkW;
dat.short_isi = short_isi; 
dat.lv = lv;
dat.cv = cv;
dat.id = data.spikesID;
dat.avg_fr = avg_fr;
dat.asymmetry = asymmetry;
end

function waves = get_waveform(temp)
for w = 1:length(temp.avgwave)
    temp.avgwave{w} = interp1(linspace(1,64,length(temp.avgwave{w})),...
        vertcat(temp.avgwave{w})',1:64)';
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


function spkW = get_spike_wavelet(waveforms,peakToTrough)

for a = 1:size(waveforms,1)
    w = waveforms(a,:)';
    if isnan(peakToTrough(a))
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
df.flipped = dat.flipped;
df.peak_to_trough = dat.peakToTrough;
df.trough_to_peak = dat.troughToPeak;
df.spk_wavelet = dat.spkW;
df.local_variation = dat.lv';
df.coeff_variation = dat.cv'; 
df.short_isi = dat.short_isi'; 
df.avg_fr = dat.avg_fr';
df.asymmetry = dat.asymmetry';
% add waveforms
df.waves = dat.waves;
end
