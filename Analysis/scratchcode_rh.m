%scratchcode_rh
i=6;
ns=2;

colorcode='HD';
tetrode=strsplit(data.spikesID.paths{i},filesep);
tetrode=tetrode{end};
trodeID=str2double(extractBetween(tetrode,'TT','.'));

[data_video_spk,~]=createframes_w_spikebinary(data,ns,i);
figure;
ax=gca;
subplot(1,3,1)
postprocessFigures.plot_HD_tuning(data,ns,i)
subplot(1,3,2)
postprocessFigures.spikesonpath_2d(ax,data_video_spk,data.lfp.ts,data.lfp.theta_phase(trodeID,:),colorcode)
% PLOT RATEMAP
subplot(1,3,3)
ax=gca
postprocessFigures.ratemaps_2d(ax,data.ratemap{i,ns},data.measures(i,:,ns),data.varnames)
dir_bins=0:360/8:360;
for i=1:length(dir_bins)-1
    dir_bin_dat{i}=data_video_spk(data_video_spk(:,4)>=dir_bins(i) & data_video_spk(:,4)<dir_bins(i+1),:); 
end
figure
ax=gca;
for i=1:length(dir_bin_dat)
    subplot(4,4,i)
    postprocessFigures.spikesonpath_2d(ax,dir_bin_dat{i},data.lfp.ts,data.lfp.theta_phase(trodeID,:),colorcode)
end




theta=0:.01:2*pi;
color=hsv(length(theta));
figure;
subplot(2,1,1)
plot(data_video_spk(:,1),data_video_spk(:,2),'.k')
hold on
spkbinary=logical(data_video_spk(:,6));
scatter(data_video_spk(spkbinary,1),data_video_spk(spkbinary,2),20,...
    interp1(rad2deg(theta)',color,data_video_spk(spkbinary,4)),'filled');
ylabel('x axis')

subplot(2,1,2)
plot(data_video_spk(:,1),data_video_spk(:,3),'.k')
hold on
spkbinary=logical(data_video_spk(:,6));
scatter(data_video_spk(spkbinary,1),data_video_spk(spkbinary,3),20,...
    interp1(rad2deg(theta)',color,data_video_spk(spkbinary,4)),'filled');
ylabel('y axis')
xlabel('time')
darkBackground(gcf,[0.2 0.2 0.2],[0.7 0.7 0.7])


figure
scatter3(data_video_spk(:,1),data_video_spk(:,2),data_video_spk(:,3),'.k')
hold on;
spkbinary=logical(data_video_spk(:,6));
scatter3(data_video_spk(spkbinary,1),data_video_spk(spkbinary,2),data_video_spk(spkbinary,3),20,...
    interp1(rad2deg(theta)',color,data_video_spk(spkbinary,4)),'filled');
xlabel('time')
ylabel('x axis')
zlabel('y axis')
darkBackground(gcf,[0.2 0.2 0.2],[0.7 0.7 0.7])



time_bins=data_video_spk(1,1):120:data_video_spk(end,1);
clear time_bin_dat
for i=1:length(time_bins)-1
    time_bin_dat{i}=data_video_spk(data_video_spk(:,1)>=time_bins(i) & data_video_spk(:,1)<time_bins(i+1),:);
     
end

figure
for i=1:length(time_bin_dat)
    subplot(3,3,i)
    postprocessFigures.spikesonpath_2d(ax,time_bin_dat{i},data.lfp.ts,data.lfp.theta_phase(trodeID,:),colorcode)
    title([num2str(time_bins(i)),' to ',num2str(time_bins(i+1)),' sec'])
end
darkBackground(gcf,[0.2 0.2 0.2],[0.7 0.7 0.7])

figure
for i=1:length(time_bin_dat)
    subplot(3,3,i)
    
    [ ratemap,~,~,~,~] =...
        bindata(time_bin_dat{i}(time_bin_dat{i}(:,6)==0,:),30,...
        time_bin_dat{i}(time_bin_dat{i}(:,6)==1,:),0,100);
    imAlpha=ones(size(ratemap));
    imAlpha(isnan(ratemap))=0;
    imagesc(ratemap,'AlphaData',imAlpha);
    axis xy; axis off; hold on; box off; axis image;
    colormap(gca,viridis(255))
    title([num2str(time_bins(i)),' to ',num2str(time_bins(i+1)),' sec'])
end





%% concatAnalyzeResults to single table
% run after SPLIT BY REGION section in AnalyzeResults_HPCatn

data=[group1ca1;group2ca1;group1ca3;group2ca3;group1cortex;group2cortex]; 




group1ca1id = [group1ca1id,cellstr(num2str(group1ca1(:,end)))];
group2ca1id = [group2ca1id,cellstr(num2str(group2ca1(:,end)))];
group1ca3id = [group1ca3id,cellstr(num2str(group1ca3(:,end)))];
group2ca3id = [group2ca3id,cellstr(num2str(group2ca3(:,end)))];
group1cortexid = [group1cortexid,cellstr(num2str(group1cortex(:,end)))];
group2cortexid = [group2cortexid,cellstr(num2str(group2cortex(:,end)))];

tempid=[group1ca1id;group2ca1id;group1ca3id;group2ca3id;group1cortexid;group2cortexid];
    
% region
data(ismember(arrayfun(@(i1)fullfile((tempid{i1,:})),(1:size(tempid,1))','un',0),...
    arrayfun(@(i1)fullfile((group1ca1id{i1,:})),(1:size(group1ca1id,1))','un',0)),end+1)=1;

data(ismember(arrayfun(@(i1)fullfile((tempid{i1,:})),(1:size(tempid,1))','un',0),...
    arrayfun(@(i1)fullfile((group2ca1id{i1,:})),(1:size(group2ca1id,1))','un',0)),end)=1;

data(ismember(arrayfun(@(i1)fullfile((tempid{i1,:})),(1:size(tempid,1))','un',0),...
    arrayfun(@(i1)fullfile((group1ca3id{i1,:})),(1:size(group1ca3id,1))','un',0)),end)=2;

data(ismember(arrayfun(@(i1)fullfile((tempid{i1,:})),(1:size(tempid,1))','un',0),...
    arrayfun(@(i1)fullfile((group2ca3id{i1,:})),(1:size(group2ca3id,1))','un',0)),end)=2;

data(ismember(arrayfun(@(i1)fullfile((tempid{i1,:})),(1:size(tempid,1))','un',0),...
    arrayfun(@(i1)fullfile((group1cortexid{i1,:})),(1:size(group1cortexid,1))','un',0)),end)=3;

data(ismember(arrayfun(@(i1)fullfile((tempid{i1,:})),(1:size(tempid,1))','un',0),...
    arrayfun(@(i1)fullfile((group2cortexid{i1,:})),(1:size(group2cortexid,1))','un',0)),end)=3;

varnames=[varnames,'brain_region'];

% group
data(ismember(arrayfun(@(i1)fullfile((tempid{i1,:})),(1:size(tempid,1))','un',0),...
    arrayfun(@(i1)fullfile((group1ca1id{i1,:})),(1:size(group1ca1id,1))','un',0)),end+1)=1;

data(ismember(arrayfun(@(i1)fullfile((tempid{i1,:})),(1:size(tempid,1))','un',0),...
    arrayfun(@(i1)fullfile((group2ca1id{i1,:})),(1:size(group2ca1id,1))','un',0)),end)=2;

data(ismember(arrayfun(@(i1)fullfile((tempid{i1,:})),(1:size(tempid,1))','un',0),...
    arrayfun(@(i1)fullfile((group1ca3id{i1,:})),(1:size(group1ca3id,1))','un',0)),end)=1;

data(ismember(arrayfun(@(i1)fullfile((tempid{i1,:})),(1:size(tempid,1))','un',0),...
    arrayfun(@(i1)fullfile((group2ca3id{i1,:})),(1:size(group2ca3id,1))','un',0)),end)=2;

data(ismember(arrayfun(@(i1)fullfile((tempid{i1,:})),(1:size(tempid,1))','un',0),...
    arrayfun(@(i1)fullfile((group1cortexid{i1,:})),(1:size(group1cortexid,1))','un',0)),end)=1;

data(ismember(arrayfun(@(i1)fullfile((tempid{i1,:})),(1:size(tempid,1))','un',0),...
    arrayfun(@(i1)fullfile((group2cortexid{i1,:})),(1:size(group2cortexid,1))','un',0)),end)=2;

varnames=[varnames,'group_id'];

tempid(:,end)=[];

varnames=['rat_session','tt','cell',varnames];

varnames = regexprep(varnames, '\s', '');


T = array2table(horzcat(tempid,num2cell(data)));
T.Properties.VariableNames = varnames;

writetable(T,...
   'D:\Dropbox\school work\UNM\Classes\Stats_527\HPCatn_data.csv'); %save data


%% delete .dat files
cd F:\ryan_harvey\Projects\PAE_PlaceCell\data
files=dir('**/*.dat');
for i=1:size(files,1)
    disp(fullfile(files(i).folder,files(i).name))
    delete(fullfile(files(i).folder,files(i).name));
end
disp([num2str(sum([files(:).bytes])/1000000000), 'GBs deleted'])


%% Rasters
data = load('F:\Projects\ATNppc\ProcessedData\ATNppc01_S20190909163020.mat') 
cortex=[1,2,3,4,5,14,15,16];
thalamus=[6,7,8,9,10,11,12,13];

tetrodes = data.spikesID.TetrodeNum ;
tetrodes = extractBetween(tetrodes,'TT','.mat');

cortexidx = find(ismember(str2double([tetrodes]),cortex));
thalamusidx = find(ismember(str2double([tetrodes]),thalamus));

X = cellfun(@transpose,data.Spikes,'un',0);
% X = X([cortexidx;thalamusidx]);

figure
subplot(2,1,1)
LineFormat = struct()
LineFormat.Color = [0.9100    0.4100    0.1700];
cortex_spikes = X(cortexidx);
[~,I] = sort(cellfun(@length,cortex_spikes));
cortex_spikes = cortex_spikes(I);
[x1,y1]=plotSpikeRaster(cortex_spikes,'PlotType','vertline','LineFormat',LineFormat);hold on;

subplot(2,1,2)
LineFormat = struct()
LineFormat.Color = [.2 .2 .2];
thalamus_spikes = X(thalamusidx);
[~,I] = sort(cellfun(@length,thalamus_spikes));
thalamus_spikes = thalamus_spikes(I);
[x2,y2]=plotSpikeRaster(thalamus_spikes,'PlotType','vertline','LineFormat',LineFormat);hold on;
[0.9100    0.4100    0.1700]

fig = figure; fig.Color = [1 1 1];
plot(x1,y1+max(y2),'Color',[0.9100    0.4100    0.1700]);hold on
plot(x2,y2,'Color',[0 0 1])
% darkBackground(gcf,[0.1 0.1 0.1],[0.7 0.7 0.7])
axis tight
xlabel('Time (sec)')
ylabel('Cells')
export_fig('F:\Projects\ATNppc\Figures\colorraster_S20190909163020.png','-m4') 


X = cellfun(@transpose,data.Spikes,'un',0);
[~,I] = sort(cellfun(@length,X));
X = X(I);
clear x y
for i=1:length(X)
    [x{i},y{i}]=plotSpikeRaster(X(i),'PlotType','vertline');
end
for i=1:length(x)
    plot(x{i},y{i}+i,'Color',[rand(1),rand(1),rand(1)]);hold on
end
axis tight
grid on
xlabel('Time (sec)')
ylabel('Cells')


%% bootstrap std rate maps

mins = 6;
n = mins*60;
% bins = data_video_nospk(1,1):120:data_video_nospk(end,1);
for ibin = 1:1000
    %     idx = data_video_spk(:,1) >= bins(ibin) & data_video_spk(:,1) <= bins(ibin+1);
    %     idx = randi([0 1],length(data_video_spk),1);
    
    idx = zeros(length(data_video_spk),1);
    idx(randperm(numel(idx), n)) = 1;
    
    [ratemap(:,:,ibin),nBinsx,nBinsy,occ,Coherence]=bindata(data_video_spk(idx & data_video_spk(:,6)==0 ,:),...
        data.samplerate,...
        data_video_spk(idx & data_video_spk(:,6)==1,:),...
        track,...
        data.maze_size_cm(event));
end

figure
subplot(1,3,1)
pcolor(SmoothRateMap2)
axis image
hold on
title('rate map')
colorbar 

subplot(1,3,2)
pcolor(nanmean(ratemap,3))
axis image
hold on
title('mean bootstrapped rate map')
colorbar 


subplot(1,3,3)
pcolor(nanstd(ratemap,[],3))
axis image
hold on
title('std')
colorbar 

%% draw max diameter
D = pdist([x_temp',y_temp']);
[fieldwidth,I] = max(D);

Z = squareform(D);
[row,col] = find(Z==fieldwidth);


Z = squareform(D);
[row2,col2] = find(Z==D(I+1));

figure
plot(x_temp([row2(1),row(1)]),y_temp([row2(2),row(2)]))
hold on
plot(x_temp,y_temp)