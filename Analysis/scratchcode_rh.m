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


%% use bandpass filtered waveforms 
            
myKsDir='F:\Projects\PAE_PlaceCell\data\LEM3216\2019-08-23_13-56-19';
sp.n_channels_dat=64;

datfile=strsplit(myKsDir,filesep);datfile=[datfile{end},'.dat'];
            gwfparams.dataDir = [myKsDir,filesep];    % KiloSort/Phy output folder
            gwfparams.fileName = datfile;         % .dat file containing the raw
            gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
            gwfparams.nCh = sp.n_channels_dat;        % Number of channels that were streamed to disk in .dat file

% Load .dat and KiloSort/Phy output
fileName = fullfile(gwfparams.dataDir,gwfparams.fileName);   
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh nSamp], 'x'});

chMap = readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy'))+1;% Order in which data was streamed to disk; must be 1-indexed for Matlab
nChInMap = numel(chMap);

fid = fopen('filtered.dat', 'w');

% [b1,a1] = butter(3,[500/(32000/2) (32000*.475)/(32000/2)])
[b1, a1] = butter(3, 150/32000*2, 'high')

for i=1:nChInMap
    
    dataRAW = mmf.Data.x(i,:);
            
    dataRAW = filter(b1, a1, dataRAW);
    dataRAW = flipud(dataRAW);
    dataRAW = filter(b1, a1, dataRAW);
    dataRAW = flipud(dataRAW);
    
    fwrite(fid, dataRAW, 'int16');
    
end
fclose(fid);


for i=1:64
    dataRAW = mmf.Data.x(25,:);
    
    
    datr = filter(b1, a1, dataRAW);
    datr = flipud(datr);
    datr = filter(b1, a1, datr);
    datr = flipud(datr);
    figure;plot(datr(1000:2000))
    pause(.0001)
end


fileName = 'filtered.dat';   
gwfparams.nCh=64;
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
mmf_filt = memmapfile(fileName, 'Format', {gwfparams.dataType, [nSamp gwfparams.nCh  ], 'x'});

figure;plot( mmf_filt.Data.x(1000:2000,25))




fid =fopen('test.dat','w')
for i=1:10
        fwrite(fid, [1:100]*i, 'int16');
 
end
fclose(fid);

fileName = 'test.dat';   
gwfparams.nCh=100;
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
mmf_test = memmapfile(fileName, 'Format', {'int16', [gwfparams.nCh nSamp], 'x'});






myKsDir='F:\Projects\HPCatn\data\HPCatn07\2019-09-09_12-22-06';
sp.n_channels_dat=64;

datfile=strsplit(myKsDir,filesep);datfile=[datfile{end},'.dat'];
            gwfparams.dataDir = [myKsDir,filesep];    % KiloSort/Phy output folder
            gwfparams.fileName = datfile;         % .dat file containing the raw
            gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
            gwfparams.nCh = sp.n_channels_dat;        % Number of channels that were streamed to disk in .dat file

fileName = fullfile(gwfparams.dataDir,gwfparams.fileName);   
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel

% open raw file
fid_raw = fopen('2019-09-09_12-22-06.dat','r');
% open filtered file to write to
fid_filt = fopen('filtered_2.dat', 'w'); % open .dat file for writing
ntbuff =65600;
Nbatch      = ceil(nSamp/(ntbuff));
[b1, a1] = butter(3, 150/32000*2, 'high');

for ibatch = linspace(0,nSamp,Nbatch)
    fprintf('batch %3.0f of %3.0f \n', ibatch,nSamp);

    fseek(fid_raw, ibatch, 'bof');
    % read raw file
    buff=fread(fid_raw,[64,ntbuff],'*int16');
    
    dataRAW = gpuArray(buff);
    dataRAW = dataRAW';
    dataRAW = single(dataRAW);
    
    % subtract the mean from each channel
    dataRAW = dataRAW - mean(dataRAW, 1);    
    
    datr = filter(b1, a1, dataRAW);
    datr = flipud(datr);
    datr = filter(b1, a1, datr);
    datr = flipud(datr);
    
    datr = datr - median(datr, 2);
    
    % write file
    datcpu  = gather_try(int16(datr));
    fwrite(fid_filt, datcpu, 'int16');
 end
fclose(fid_raw);
fclose(fid_filt);
