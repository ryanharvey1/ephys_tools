% Compile Attract compiles results from attractor_xcor scripts. 
% LB & RH 10/2020
%   
%
%% 1.) Import HD cell list and velocity results (includes cell_idx/cell types). 
data_path = 'F:\ClarkP30_Recordings\Analysis\Attractor\';
hd_list = readtable('F:\ClarkP30_Recordings\Analysis\hd_cell_list.csv');
hd_cells = readtable('F:\ClarkP30_Recordings\Analysis\hd_cell_spike_idx.csv');
hd_data = readtable('F:\ClarkP30_Recordings\Analysis\hd_data_all.csv');
vel_class = readtable('F:\ClarkP30_Recordings\Analysis\Velocity\velocity_res.csv');

% Remove remove conditions beyond baseline for hd_data_all (will retrieve preferred direction from here). 
hd_data = hd_data(hd_data.Condition_num == 1,:);

%% 2.) Build pairwise data frame to store cell-pair information

% Retrieve attractor_cor data structures 
files = dir([data_path ,'*.mat']);
% initalize matrices 
df = []; name = []; xCor =[]; space_cor = []; central_sec = [];

% Loop to grab data from attractor_res
session = 1;
for i = 1:length(files)
   % load attractor_res file
   load([data_path,files(i).name]);
   hd_idx = zeros(size(res{1,1}.pair_ID,1),1);
      
   if length(res)  ~= 4
       continue
   end
   
   % Concatenate IDs, r, and central 1s for first session
   df = [df; res{1,session}.pair_ID, res{1,session}.r', res{1,session}.central_sec'];
   name = [name;  repmat(extractAfter(files(i).name,'attractor_res_'),size(res{1,session}.r,2),1)];
   xCor = [xCor; res{1,session}.zcor];
   space_cor = [space_cor; res{1,session}.spatial_corr];
   central_sec = [central_sec; res{1,session}.central_sec'];
   clear res
end

% Normalize space_cor to make comparison easier
space_cor = space_cor./sum(space_cor,2);
lags = -20:.2:20; 

% Save Dataframe as table 
to_save = table(name,df(:,1),df(:,2),df(:,3),num2cell(central_sec),'VariableNames',...
    {'SessionID','ref_cell','test_cell','r','central_sec'});

writetable(to_save,[data_path ,'attractor_estimates.csv']);
clear to_save 

%% 3.) Reload pair-wise data frame aka Attractor_estimates
attractor_estimates = readtable([data_path ,'attractor_estimates.csv']);

attractor_estimates.SessionID = extractBefore(attractor_estimates.SessionID,'.mat');
hd_list.SessionID = extractBefore(hd_list.SessionID,'.mat');

%% 4.) Add cell types to Attractor_estimates table 
for i = 1:length(attractor_estimates.SessionID)
    % loop session index
    cur_ses = attractor_estimates.SessionID{i};
    ref_idx = attractor_estimates.ref_cell(i);
    test_idx = attractor_estimates.test_cell(i);
    ref_cell_type(i,1) = vel_class.cell_type(ismember(vel_class.SessionID,cur_ses) & vel_class.cell_idx == ref_idx,1);
    test_cell_type(i,1) = vel_class.cell_type(ismember(vel_class.SessionID,cur_ses) & vel_class.cell_idx == test_idx,1);
end

% Add to table here
attractor_estimates.ref_type = ref_cell_type;
attractor_estimates.test_type = test_cell_type;

%% 5.) Get index for all HD cell pairs 
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

big_idx = logical(big_idx);
group = contains(attractor_estimates.SessionID,'LB03'); % add grouping variable

% Create index for all int-int pyr-int/int-pyr, and pyr-pyr pairs; 
attractor_estimates.int_pair = ismember(attractor_estimates.ref_type,'int') & ismember(attractor_estimates.test_type,'int');
attractor_estimates.pyr_pair = ismember(attractor_estimates.ref_type,'pyr') & ismember(attractor_estimates.test_type,'pyr');
attractor_estimates.pyr_int_pair = ismember(attractor_estimates.ref_type,'pyr') & ismember(attractor_estimates.test_type,'int');
attractor_estimates.int_pyr_pair = ismember(attractor_estimates.ref_type,'int') & ismember(attractor_estimates.test_type,'pyr');

%% 6.) Create Pairs table of just cells that met for HD 

% Start by choosing only pairs where both met for HD 
pairs_table = attractor_estimates(big_idx,:);

% Keep HD pairs from xCOR and Space cor
xCor_cor_hd = xCor(big_idx,:);
space_cor_hd = space_cor(big_idx,:);
group_hd = group(big_idx,:);

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

% Add angular difference and cell types to cell pairs table
pairs_table.ang_diff = wrapTo180(rad2deg(circ_dist(deg2rad(pd_ref),deg2rad(pd_test))));

%% 7.) Check rayleigh test for HD cell pairs 

% keep cells that meet for r-test with p > .001 and Z > 5 
[ref_pval, ref_z, ref_group,ref_tuning] = batch_rtest(pairs_table.SessionID,pairs_table.ref_cell);
[test_pval, test_z, test_group, test_tuning] = batch_rtest(pairs_table.SessionID,pairs_table.test_cell);

keep_idx = ref_pval < .001 & test_pval < .001 & ref_z > 5 & test_z > 5;

% Index tuning curves, xcor, space cor and group index for highly directional HD cells. 
keepHD = test_tuning(keep_idx,:);
xCor_cor_hd = xCor_cor_hd(keep_idx,:);
space_cor_hd = space_cor_hd(keep_idx,:);
group_hd = group_hd(keep_idx,:);
pairs_table = pairs_table(keep_idx,:);

%% PLOTS BELOW

%% Non-HD versus highly HD pairs 

summary_plot(attractor_estimates,pairs_table,~big_idx,group_hd,'F344','TgF344-AD',{'Non-HD pairs','HD pairs'})

%% Pyramidal-Interneuron and Interneuron-Pyramidal Pairs 


%% Interneuron-Interneuron vs Pyramidal-Pyramidal Pairs Summary 


%Grab xCor for each animal
wt_xCor = xCor(big_idx & group == 1,:);
tg_xCor = xCor(big_idx & group == 0,:);

wt_space = space_cor_hd(big_idx & group == 1,:);
tg_space = space_cor_hd(big_idx & group == 0,:);

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

function smoothed = smooth_hd(tuning)
gaus = @(tuning,mu,sig)exp(-(((tuning-mu).^2)/(2*sig.^2)));
window = 20;
x = -window/2:window/2;
mu = 0;
sig = 3;
kernel = gaus(x,mu,sig);
smoothed = conv(tuning,kernel/sum(kernel),'same');
end

function summary_plot(pairs_1,pairs_2,group_idx_1,group_idx_2,group_1,group_2,pairs_name)

% pairs_1: table for cell pairs of interest with r & central_sec variables
% pairs_2: table for cell pairs of interest with r & central_sec variables
% group_idx: logical index of size of pairs_1 or pairs_2 respectively
% group_1: cell array of group1 name i.e. 'F344'
% group_2: cell array of group1 name i.e. 'TgF344'
% pairs_name: cell array for legened {'HD:HD','Non-HD:HD'}

edges = 0:.05:max([pairs_1.r;pairs_2.r]);
fig = figure; 
fig.Color = [1 1 1];
subplot(2,2,1)
histogram(pairs_1.r(group_idx_1 == 1),80,'Normalization','prob','BinEdges',edges)
hold on;
histogram(pairs_2.r(group_idx_2 == 1),80,'EdgeAlpha',.5, 'EdgeColor','r',...
'FaceAlpha',.5, 'FaceColor','r','Normalization','prob','BinEdges',edges)
xlabel('Rayleigh Vector')
title(['R Length of Spatial xCor:',group_1])
ylabel('Probability')
legend(pairs_name)
set(gca,'FontSize',12,'FontWeight','bold')
subplot(2,2,2)
histogram(pairs_1.r(group_idx_1 == 0),80,'Normalization','prob','BinEdges',edges)
hold on;
histogram(pairs_2.r(group_idx_2 == 0),80,'EdgeAlpha',.5, 'EdgeColor','r',...
'FaceAlpha',.5, 'FaceColor','r','Normalization','prob','BinEdges',edges)
set(gca,'FontSize',12,'FontWeight','bold')
title(['R Length of Spatial xCor:',group_2])

edges = min([pairs_1.central_sec;pairs_2.central_sec]):.05:max([pairs_1.central_sec;pairs_2.central_sec]);
subplot(2,2,3)
histogram(pairs_1.central_sec(group_idx_1 == 1),100,'Normalization','prob','BinEdges',edges)
hold on;
histogram(pairs_2.central_sec(group_idx_2 == 1),100,'EdgeAlpha',.5, 'EdgeColor','r',...
'FaceAlpha',.5, 'FaceColor','r','Normalization','prob','BinEdges',edges)
xlabel('Normalized Correlation')
ylabel('Probability')
legend(pairs_name)
title(['Mean Central Second xCor:', group_1])
set(gca,'FontSize',12,'FontWeight','bold')
subplot(2,2,4)
histogram(pairs_1.central_sec(group_idx_1 == 0),100,'Normalization','prob','BinEdges',edges)
hold on;
histogram(pairs_2.central_sec(group_idx_2 == 0),100,'EdgeAlpha',.5, 'EdgeColor','r',...
'FaceAlpha',.5, 'FaceColor','r','Normalization','prob','BinEdges',edges)
set(gca,'FontSize',12,'FontWeight','bold')
title(['Mean Central Second xCor:', group_2])
end
