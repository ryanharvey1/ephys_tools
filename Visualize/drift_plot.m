% set figure defaults

function drift_plot(drift_files)
fig=figure('DefaultAxesFontSize',8,'defaultAxesFontName','sans-serif','defaultTextFontName','sans-serif');
fig.Color = [1,1,1];
[fig_width_in, fig_height_in] = set_size('thesis', 1.5, [3,4]);
set(fig,'Position',[1 3 fig_width_in fig_height_in])
SessionID = drift_files.SessionID{1};

data = load(['F:\ClarkP30_Recordings\ProcessedData\',drift_files.SessionID{1}],'spikesID');
cell = drift_files.cell(1);
tet = sscanf(drift_files.tetrode{1},'TT%d.mat');
cell_idx = find_cells(data,tet,cell);


conditions = {'Baseline','Standard 1','Rotation','Standard 2'};
% plot tuning curves
subplot(3,4,1)
example_HD(SessionID,1,cell_idx)
title(conditions(1),'FontWeight','normal')
subplot(3,4,2)
example_HD(SessionID,2,cell_idx)
title(conditions(2),'FontWeight','normal')
subplot(3,4,3)
example_HD(SessionID,3,cell_idx)
title(conditions(3),'FontWeight','normal')
subplot(3,4,4)
example_HD(SessionID,4,cell_idx)
title(conditions(4),'FontWeight','normal')
% plot moving turning curves over time bins
subplot(3,4,5)
moving_tuning(SessionID,1,cell_idx,1)
subplot(3,4,6)
moving_tuning(SessionID,2,cell_idx,0)
subplot(3,4,7)
moving_tuning(SessionID,3,cell_idx,0)
subplot(3,4,8)
moving_tuning(SessionID,4,cell_idx,0)
subplot(3,4,9)
hd_tune_over_time(SessionID,1,cell_idx,1)
subplot(3,4,10)
hd_tune_over_time(SessionID,2,cell_idx,0)
subplot(3,4,11)
hd_tune_over_time(SessionID,3,cell_idx,0)
subplot(3,4,12)
hd_tune_over_time(SessionID,4,cell_idx,0)

end

function moving_tuning(SessionID,cond,cell,set_ylabel)
data = load(['F:\ClarkP30_Recordings\ProcessedData\',SessionID],'frames','Spikes','spikesID','events');
idx = data.frames(:,1) >= data.events(1,cond) & data.frames(:,1) <= data.events(2,cond);
data.frames = data.frames(idx,:);
ts_bins = min(data.frames(:,1)):20:max(data.frames(:,1));
for i = 1:length(ts_bins)-1
    temp_frames = data.frames(data.frames(:,1)>=ts_bins(i) & data.frames(:,1)<=ts_bins(i+1),:);
    ts = temp_frames(:,1);
    ang = temp_frames(:,4);
%     tet = sscanf(data.spikesID.TetrodeNum{cell},'TT%d.mat');
%     cell_idx = find_cells(data,tet,data.spikesID.CellNum(cell));
    spikes = data.Spikes{cell,1};
    spikes = spikes(spikes>=ts_bins(i) & spikes<=ts_bins(i+1));
    % Get angle spike
    try
    angspk = circular_interp(ts,ang,spikes);
    catch
        smoothed(i,:) = zeros(size(hdTuning,2),1);
        continue
    end
    [r,I,Ispk,peakrate,prefdirec,hdTuning]=tuningcurve(ang,angspk,30);
    smoothed(i,:) = smooth_hd(hdTuning);
end
imagesc(smoothed')
colormap(magma)
axis xy
set(gca,'xtick',[])
set(gca,'ytick',[20,40,60],'yticklabels',[120,240,360])
if set_ylabel
    ylabel('degrees')
else
    axis off
end
% hold on;
% [m,y] = max(smoothed');
% x = 1:size(smoothed,1);
% idx = m == 0;
% 
% x(idx) = [];
% y(idx) = [];
% 
% % find(abs(diff(y))*6 > 180)
% 
% plot(x,y,'w','LineWidth',.5)

end
% plot head direction over time 
function hd_tune_over_time(SessionID,cond,cell,set_label)
data = load(['F:\ClarkP30_Recordings\ProcessedData\',SessionID],'frames','Spikes','spikesID','events');
% filter by cond
idx = data.frames(:,1) >= data.events(1,cond) & data.frames(:,1) <= data.events(2,cond);
data.frames = data.frames(idx,:);
ts = data.frames(:,1);
ang = data.frames(:,4);
% tet = sscanf(data.spikesID.TetrodeNum{cell},'TT%d.mat');
% cell_idx = find_cells(data,tet,data.spikesID.CellNum(cell));
spikes = data.Spikes{cell,1};
% Get angle spike
angspk = circular_interp(ts,ang,spikes);
% set up colormap
theta=0:.01:2*pi;
color = hsv(length(theta));
theta = theta';
plot(ts,ang,'Color',[.8 .8 .8]);
hold on
h=scatter(spikes,angspk,5,interp1(rad2deg(theta),color,angspk,'nearest'),'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
% set(h,'FaceAlpha',.5);
box off
% xlabel('Time (sec)')
% ylabel('degrees')
% axis tight
ylim([0,360])
xlim([min(ts),max(ts)])
set(gca,'ytick',[120,240,360],'yticklabels',[120,240,360])
if set_label
    ylabel('degrees')
    xlabel('Time (sec)')
else
    set(gca,'ytick',[])
end
end
function example_HD(SessionID,cond,cell)
data = load(['F:\ClarkP30_Recordings\ProcessedData\',SessionID],'hdTuning','spikesID');
% tet = sscanf(data.spikesID.TetrodeNum{cell},'TT%d.mat');
% cell_idx = find_cells(data,tet,data.spikesID.CellNum(cell));
tuning = data.hdTuning{cell,cond};
smoothed = smooth_hd(tuning);
% r = data.measures(cell,contains(data.varnames,'mean_vector_length'),cond);
% Ispk = data.measures(cell,contains(data.varnames,'Direct_infoContent'),cond);
theta = 0:.01:2*pi;
color = hsv(length(theta));
% fig=figure('Name',[data.rat,'  ',data.sessionID,'  ', 'Cell: ',num2str(data.spikesID.CellNum(cell))],'NumberTitle','off');
% fig.Color = [1 1 1];
% subaxis(2,1,1);
angBins = 0:6:360;
hd_centers = movmedian(angBins,2);
hd_centers(1)=[];
Polarplot = polar(deg2rad(hd_centers),smoothed,'b');
set(Polarplot,'linewidth',1,'color','k');
axis off
set(0,'Showhiddenhandles','on')
extrastuff = setdiff(get(gca,'children'),Polarplot);
delete(extrastuff)
hold on
horizontal=line([-max(smoothed) max(smoothed)],[0 0]); % for running max and min
vertical=line([0 0],[-max(smoothed) max(smoothed)]);
set(horizontal,'linewidth',2,'color',[.4 .4 .4]);
set(vertical,'linewidth',2,'color',[.4 .4 .4]);
axis image
uistack(horizontal,'bottom')
uistack(vertical,'bottom')
% title(sprintf('r: %4.2f DIC: %4.2f' ,[r,Ispk]))
theta=0:.01:2*pi;
color=hsv(length(theta));
[~,I] = max(smoothed);
h=fill(get(Polarplot, 'XData'), get(Polarplot, 'YData'),...
    interp1(rad2deg(theta)',color,hd_centers(I)),'nearest');
set(h,'FaceAlpha',.5)
% subaxis(2,1,2);
% postprocessFigures.avg_waveforms(data,cond,cell)
% max_lag = 0.05;
% t_bin=0.002;
% % Acor - taken from intrinsic frequency 2
% if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')
%     max_lag = t_bin*floor(max_lag/t_bin)+max_lag*t_bin;
% end
% ses_idx = data.Spikes{cell,1} > data.events(1,cond) & data.Spikes{cell,1} < data.events(2,cond);
% [cor, lag] = CrossCorr(data.Spikes{cell,1}(ses_idx),data.Spikes{cell,1}(ses_idx),'lag',...
%     [-max_lag max_lag],'binsize',t_bin,'norm','prob');
% bar(lag,cor,'k')
% xlabel('Time (ms)')
% box off
% set(gcf,'OuterPosition',[1856 361 338 605])
end
function smoothed = smooth_hd(tuning)
tuning = [tuning tuning tuning];
gaus = @(tuning,mu,sig)exp(-(((tuning-mu).^2)/(2*sig.^2)));
window = 10;
x = -window/2:window/2;
mu = 0;
sig = 3;
kernel = gaus(x,mu,sig);
smoothed = conv(tuning,kernel/sum(kernel),'same');
smoothed = smoothed(:,61:120);
end