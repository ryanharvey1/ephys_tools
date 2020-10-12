% Compile Attract compiles results from attractor_xcor scripts. 
% LB & RH 10/2020
% To-do: 
%   - Apply Rayleigh test as inclusion for attractor analyses
%   
%
%% import the files
data_path = 'F:\ClarkP30_Recordings\Analysis\Attractor\';
hd_list = readtable('F:\ClarkP30_Recordings\Analysis\hd_cell_list.csv');
hd_cells = readtable('F:\ClarkP30_Recordings\Analysis\hd_cell_spike_idx.csv');
hd_data = readtable('F:\ClarkP30_Recordings\Analysis\hd_data_all.csv');
hd_data = hd_data(hd_data.Condition_num == 1,:);

%% Build pairwise data frame
files = dir([data_path ,'*.mat']);
df = []; name = []; xCor =[]; space_cor = []; 

for i = 1:length(files)
   % load attractor_res file
   load([data_path,files(i).name]);
   hd_idx = zeros(size(res.pair_ID,1),1);
      
   % Concatenate IDs, r, and central 1s 
   df = [df; res.pair_ID, res.r', res.central_sec'];
   name = [name;  repmat(extractAfter(files(i).name,'attractor_res_'),size(res.r,2),1)];
   xCor = [xCor; res.cor];
   space_cor = [space_cor; res.spatial_corr];
    
end

%% Normalize 
space_cor = space_cor./sum(space_cor,2);
xCor = xCor./sum(xCor,2);
lags = -20:.2:20; 
central_sec = nanmean(xCor(:,lags > -0.5 & lags < 0.5),2); % normalize central second

%% Save Dataframe as table
to_save = table(name,df(:,1),df(:,2),df(:,3),num2cell(central_sec),'VariableNames',...
    {'SessionID','ref_cell','test_cell','r','central_sec'});

writetable(to_save,[data_path ,'attractor_estimates.csv']);
clear to_save 

attractor_estimates = readtable([data_path ,'attractor_estimates.csv']);

attractor_estimates.SessionID = extractBefore(attractor_estimates.SessionID,'.mat');
hd_list.SessionID = extractBefore(hd_list.SessionID,'.mat');

%% Get index for all HD cell pairs 
big_idx = zeros(length(attractor_estimates.SessionID),1);
for i = 1:length(hd_list.SessionID)
    % loop session index
    cur_ses = contains(hd_list.SessionID,hd_list.SessionID{i});
    
    % Important indicies
    ses_idx = contains(attractor_estimates.SessionID,hd_list.SessionID{i});    
    ref_idx = ismember(attractor_estimates.ref_cell,hd_cells.n(cur_ses));
    test_idx = ismember(attractor_estimates.test_cell,hd_cells.n(cur_ses));
    
    % Append indicies to list
    big_idx = big_idx + (ses_idx & ref_idx & test_idx);
end
%% Create table from all HD cell pairs 
big_idx = logical(big_idx);
group = contains(attractor_estimates.SessionID,'LB03');

pairs_table = attractor_estimates(big_idx,:);

%% Check rayleigh test for HD cell pairs 
[ref_pval, ref_z, ref_group,ref_tuning] = batch_rtest(pairs_table.SessionID,pairs_table.ref_cell);
[test_pval, test_z, test_group, test_tuning] = batch_rtest(pairs_table.SessionID,pairs_table.test_cell);

keepHD = test_tuning(ref_pval < .001 & test_pval <.001 & ref_z > 5 & test_z > 5 & ref_group == 1  & test_group == 1,:);

%%
for i = 1:length(pairs_table.SessionID)
    % loop session index
    cur_ses = pairs_table.SessionID{i};
    ref_idx = pairs_table.ref_cell(i);
    test_idx = pairs_table.test_cell(i);
    hd_idx = contains(hd_list.SessionID,pairs_table.SessionID{i});
    pd_sess = hd_data.pref_dir(hd_idx);
    pd_ref(i,1) = pd_sess(hd_cells.n(hd_idx,:) == ref_idx);
    pd_test(i,1) = pd_sess(hd_cells.n(hd_idx,:) == test_idx);
end
    
pairs_table.ang_diff = wrapTo180(rad2deg(circ_dist(deg2rad(pd_ref),deg2rad(pd_test))));
 
%% Histograms for central second and r 
fig = figure; 
fig.Color = [1 1 1];
subplot(2,2,1)
histogram(attractor_estimates.r(~big_idx& group == 1),80,'Normalization','prob')
hold on;
histogram(attractor_estimates.r(big_idx& group == 1),80,'EdgeAlpha',.5, 'EdgeColor','r',...
'FaceAlpha',.5, 'FaceColor','r','Normalization','prob')
xlabel('Rayleigh Vector')
title('R Length of Spatial xCor: F344')
ylabel('Probability')
legend({'non-hd pair','hd pair'})
set(gca,'FontSize',12,'FontWeight','bold')
subplot(2,2,2)
histogram(attractor_estimates.r(~big_idx& group == 0),80,'Normalization','prob')
hold on;
histogram(attractor_estimates.r(big_idx& group == 0),80,'EdgeAlpha',.5, 'EdgeColor','r',...
'FaceAlpha',.5, 'FaceColor','r','Normalization','prob')
set(gca,'FontSize',12,'FontWeight','bold')
title('R Length of Spatial xCor: TgF344-AD')

subplot(2,2,3)
histogram(attractor_estimates.central_sec(~big_idx& group == 1),80,'Normalization','prob')
legend({'non-hd pair','hd pair'})
hold on;
histogram(attractor_estimates.central_sec(big_idx& group == 1),80,'EdgeAlpha',.5, 'EdgeColor','r',...
'FaceAlpha',.5, 'FaceColor','r','Normalization','prob')
xlabel('Normalized Correlation')
ylabel('Probability')
title('R Length of Temporal xCor: F344')
set(gca,'FontSize',12,'FontWeight','bold')
subplot(2,2,4)
histogram(attractor_estimates.central_sec(~big_idx& group == 0),80,'Normalization','prob')
hold on;
histogram(attractor_estimates.central_sec(big_idx& group == 0),80,'EdgeAlpha',.5, 'EdgeColor','r',...
'FaceAlpha',.5, 'FaceColor','r','Normalization','prob')
set(gca,'FontSize',12,'FontWeight','bold')
title('R Length of Temporal xCor: TgF344-AD')



%Grab xCor for each animal
wt_xCor = xCor(big_idx & group == 1,:);
tg_xCor = xCor(big_idx & group == 0,:);

wt_space = space_cor(big_idx & group == 1,:);
tg_space = space_cor(big_idx & group == 0,:);

%% Spatial xCor plots 
wt_ang = pairs_table.ang_diff(contains(pairs_table.SessionID,'LB03'));
[~,wt_sort] = max(wt_ang,[],2);
[~, index] = sort(wt_ang);
wt_space_sorted  = rescale(wt_space(index, :),'InputMin',min(wt_space(index, :),[],2),'InputMax',max(wt_space(index, :),[],2));

tg_ang= pairs_table.ang_diff(~contains(pairs_table.SessionID,'LB03'));
[~,tg_sort] = max(tg_ang,[],2);
[~, index] = sort(tg_ang);
tg_space_sorted  = rescale(tg_space(index, :),'InputMin',min(tg_space(index, :),[],2),'InputMax',max(tg_space(index, :),[],2));

fig = figure; 
fig.Color = [1 1 1];
subplot(1,2,1)
imagesc(wt_space_sorted)
title('F344')
ylabel(['Cell pairs (n = ',num2str(length(wt_space_sorted)),')'])
xlabel('Angular Bins (6 degree/bin)')
set(gca,'FontSize',12,'FontWeight','bold')
subplot(1,2,2)
imagesc(tg_space_sorted)
title('TgF344-AD')
ylabel(['Cell pairs (n = ',num2str(length(tg_space_sorted)),')'])
set(gca,'FontSize',12,'FontWeight','bold')

%% Temporal xCor plots 
wt_ang = pairs_table.ang_diff(contains(pairs_table.SessionID,'LB03'));
[~, index] = sort(wt_ang);
wt_xcor_sorted  = rescale(wt_xCor(index, :),'InputMin',min(wt_xCor(index, :),[],2),'InputMax',max(wt_xCor(index, :),[],2));

tg_ang= pairs_table.ang_diff(~contains(pairs_table.SessionID,'LB03'));
[~, index] = sort(tg_ang);
tg_xcor_sorted  = rescale(tg_xCor(index, :),'InputMin',min(tg_xCor(index, :),[],2),'InputMax',max(tg_xCor(index, :),[],2));

fig = figure; 
fig.Color = [1 1 1];
subplot(1,2,1)
imagesc(wt_xcor_sorted)
title('F344')
ylabel(['Cell pairs (n = ',num2str(length(wt_xcor_sorted)),')'])
xticks(1:50:250)
xticklabels({'-20','-10','0','10','20'})
xlabel('Time lags (s)')
set(gca,'FontSize',12,'FontWeight','bold')
subplot(1,2,2)
imagesc(tg_xcor_sorted)
xticks(1:50:250)
xticklabels({'-20','-10','0','10','20'})
xlabel('Time lags (s)')
title('TgF344-AD')
ylabel(['Cell pairs (n = ',num2str(length(tg_xcor_sorted)),')'])
set(gca,'FontSize',12,'FontWeight','bold')

% xCor examples for cells with PD 0, 90, 180 
pair_cor = xCor(big_idx,:);

data = load(['F:\ClarkP30_Recordings\ProcessedData\LB10_S20200405194704.mat'],'hdTuning');
plot(data.hdTuning{7},'-k','LineWidth',2); hold on
plot(data.hdTuning{9},'-r','LineWidth',2); hold off
figure; bar(pair_cor(1118,:))

%% Local functions
function [pval, z, group,tuningcurve] = batch_rtest(files,cell)

for i = 1:length(files)
    
    data = load(['F:\ClarkP30_Recordings\ProcessedData',filesep,files{i}],'hdTuning');
    tuning = data.hdTuning{cell(i)};
    smooth_tune = smooth_hd([tuning tuning tuning]);
    smooth_tune = smooth_tune(1, 61:120);
    [pval(i,1), z(i,1)] = circ_rtest(0:2*pi/59:2*pi, smooth_tune);
    group(i,1) = ~contains(files{i},'LB03');
    tuningcurve(i,:) = smooth_tune;
end
end


%%
function smoothed = smooth_hd(tuning)
gaus = @(tuning,mu,sig)exp(-(((tuning-mu).^2)/(2*sig.^2)));
window = 20;
x = -window/2:window/2;
mu = 0;
sig = 3;
kernel = gaus(x,mu,sig);
smoothed = conv(tuning,kernel/sum(kernel),'same');
end


