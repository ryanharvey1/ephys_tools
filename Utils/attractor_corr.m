
function attractor_corr()
dataset = 'F:\ClarkP30_Recordings\ProcessedData\';
save_path = 'F:\ClarkP30_Recordings\Analysis\Attractor\';
hd = readtable('F:\ClarkP30_Recordings\Analysis\hd_cell_list.csv');

sessions = unique(hd.SessionID);
waitmessage = parfor_wait(length(sessions),'Waitbar',false,'ReportInterval',1);

for i = 1:length(sessions)
    
    % check to see if session was already run. If so, continue
    if exist([save_path,'attractor_res_',sessions{i}],'file')
        continue
    end
    
    % load data
    data = load([dataset,sessions{i}],'frames','Spikes','events');
    
    % Run Attractor analysis
    res = attractor_main(data);
    
    save_res([save_path,'attractor_res_',sessions{i}],res)
    waitmessage.Send;
end
waitmessage.Destroy

end



function save_res(save_path,res)
 save(save_path,'res')
end

function res = attractor_main(data)
% Performs spatial and temporal cross correlation analyses from Bassett,
% Wills, Cacucci 2018.
% Input:
%   data: ephys_tools data structure
%
% Output:
%   res: Structur containing 
% 
%

spikes = data.Spikes;
events = data.events(:,1);
frames = get_baseline_frames(data,events);

res = struct;
% Loop through all pairwise comparisons
cmp = [];
p = 1;
for i = 1:length(data.Spikes)
    for ii = i+1:length(data.Spikes)
        cmp = [cmp;[i, ii]];
        
        % Gather spikes
        ref = get_baseline_spikes(i,spikes,events);
        test = get_baseline_spikes(ii,spikes,events);
        
        % Calculate Spatial Correlation
        [spatial_corr(p,:), r(p)] = spatial_xcorr(ref,test,frames);
        
        % Calculate Spa
        [cor(p,:),central_sec(p)]  = temporal_xcorr(ref,test);
        p = p + 1;
    end
    
end

res.spatial_corr = spatial_corr;
res.r = r;
res.cor = cor;
res.central_sec = central_sec;
res.pair_ID = cmp;
end

%% Plot
% subplot(1,3,1)
% plot_tuning(tuning{cell_1,:})
% subplot(1,3,2)
% plot_tuning(tuning{cell_2,:})
% subplot(1,3,3)
% plot_spatial_xcorr(spatial_corr)



%% Local Functions
function ref = get_baseline_spikes(celln,spikes,events)
spike_idx = spikes{celln} > events(1) & spikes{celln} < events(2);
ref = spikes{celln}(spike_idx);

end

function [cor,central_sec] = temporal_xcorr(ref,test)

max_lag = 20;
% t_bin=0.005;
t_bin=0.200;
% Acor - taken from intrinsic frequency 2
if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')
    max_lag = t_bin*floor(max_lag/t_bin)+.5*t_bin;
end
[cor, lags] = CrossCorr(ref,test,'lag',[-max_lag max_lag],'binsize',t_bin,'norm','prob');
 cor = cor./sum(cor,2); % normalize
% Get mean 1s
central_sec = nanmean(cor(lags > -0.5 & lags < 0.5));

end

function [spatial_corr, r] = spatial_xcorr(ref,test,frames)

% subsample spikes to decrease run time 
spike_num = randperm(length(ref),min(length(ref),1000));

for iii = 1:length(spike_num)
    
    out_bound = ref(spike_num(iii)) + 10;
    in_bound = ref(spike_num(iii));
    
    index = test > in_bound & test < out_bound;
    hd_idx = frames(:,1) > in_bound & frames(:,1) < out_bound;
    current_spikes = test(index);
    angles = frames(hd_idx, 4);
    ts = frames(hd_idx,1);
    
    % Continue if ts or spikes are inadequate
    if length(ts) < 2 || isempty(current_spikes)
        continue
    end
    % line up spikes with HD
    angle_spike = interp1(ts,angles,current_spikes);
    
    % Gather occupancy and spike_angle
    [bin_angle(iii,:) , bin_spikes(iii,:)] =  binned_spikes(angles,angle_spike);
    
end

if ~exist('bin_spikes','var')
    spatial_corr = nan(1,60); 
    r = NaN;
    return
end

%% sum & smooth histograms
bin_spikes_smooth = smoothdata(sum([bin_spikes, bin_spikes, bin_spikes] ,1),2,'gaussian',5);
bin_spikes_smooth = bin_spikes_smooth(1,61:120);

bin_angle_smooth = smoothdata(sum([bin_angle, bin_angle, bin_angle] ,1),2,'gaussian',5);
bin_angle_smooth = bin_angle_smooth(1,61:120);

%% Normalize
spatial_corr = (bin_spikes_smooth./bin_angle_smooth)*30;

%% Directional Information Content of spatial cross corr
% probability of occupancy:
% Px = bin_angle_smooth./sum(bin_angle_smooth);
% logTerm = log2(spatial_corr/mean(spatial_corr));
% % Correct for undefined values
% logTerm(spatial_corr==0) = 0;
% % Little trick to express a sum as a dot product
% dic = spatial_corr * (logTerm.*Px)' ;
angBins=0:6:360;
bin_centers=movmedian(angBins,2);
bin_centers(1)=[];
r = circ_r(deg2rad(bin_centers)',spatial_corr',deg2rad(6));

end

function frames = get_baseline_frames(data,events)
idx = data.frames(:,1) > events(1) & data.frames(:,1) < events(2);
frames = data.frames(idx,:);
end

function [bin_angle , bin_spikes] =  binned_spikes(angles,angle_spike)
%angular bins
% da=pi/30;angBins=da/2:da:2*pi-da/2;
angBins = 0:6:360;
% Occupancy
bin_angle = histcounts(angles,angBins);
% Number of spikes per bin
bin_spikes = histcounts(angle_spike,angBins);
end

function plot_tuning(tuning)
smoothed_tuning = smoothdata([tuning,tuning,tuning],'gaussian',5);
angBins=0:6:360;
bin_centers=movmedian(angBins,2);
bin_centers(1)=[];
Polarplot = polar(deg2rad(bin_centers),smoothed_tuning(61:120),'b');
set(Polarplot,'linewidth',1,'color','k');
axis off
set(0,'Showhiddenhandles','on')
extrastuff = setdiff(get(gca,'children'),Polarplot);
delete(extrastuff)
hold on
horizontal=line([-max(smoothed_tuning) max(smoothed_tuning)],[0 0]); % for running max and min
vertical=line([0 0],[-max(smoothed_tuning) max(smoothed_tuning)]);
set(horizontal,'linewidth',2,'color',[.4 .4 .4]);
set(vertical,'linewidth',2,'color',[.4 .4 .4]);
axis image

uistack(horizontal,'bottom')
uistack(vertical,'bottom')

%             h=fill(get(Polarplot, 'XData'), get(Polarplot, 'YData'),...
%                 'k');
%             set(h,'FaceAlpha',.5)

% color tuning by smoothed peak direction
theta=0:.01:2*pi;
color=hsv(length(theta));
[~,I] = max(smoothed_tuning(61:120));
h=fill(get(Polarplot, 'XData'), get(Polarplot, 'YData'),...
    interp1(rad2deg(theta)',color,bin_centers(I)));
set(h,'FaceAlpha',.5)
end

function plot_spatial_xcorr(spatial_corr)
angBins=0:6:360;
bin_centers=movmedian(angBins,2);
bin_centers(1)=[];
Polarplot = polar(deg2rad(bin_centers),spatial_corr,'b');
set(Polarplot,'linewidth',1,'color','k');
axis off
set(0,'Showhiddenhandles','on')
extrastuff = setdiff(get(gca,'children'),Polarplot);
delete(extrastuff)
hold on
horizontal=line([-max(spatial_corr) max(spatial_corr)],[0 0]); % for running max and min
vertical=line([0 0],[-max(spatial_corr) max(spatial_corr)]);
set(horizontal,'linewidth',2,'color',[.4 .4 .4]);
set(vertical,'linewidth',2,'color',[.4 .4 .4]);
axis image

uistack(horizontal,'bottom')
uistack(vertical,'bottom')

%             h=fill(get(Polarplot, 'XData'), get(Polarplot, 'YData'),...
%                 'k');
%             set(h,'FaceAlpha',.5)

% color tuning by smoothed peak direction
% theta=0:.01:2*pi;
% color=hsv(length(theta));
% [~,I] = max(spatial_corr);
% h=fill(get(Polarplot, 'XData'), get(Polarplot, 'YData'),...
%     interp1(rad2deg(theta)',color,bin_centers(I)));
% set(h,'FaceAlpha',.5)
end