
%% 1. Perform shuffling procedure on processed data 

% Load cell list 
groupid = readtable('d:/Users/BClarkLab/github/Berkowitz_et_al_2021/cell_list_good.csv');
groupid = table2cell(groupid);

%% Perform direction IC shuffling on all cells
[shuff_pass,p_val,z]= shuff(groupid,'feature',{'dic'},'nshuffle',1000);
 
%% Calculate modulation depth

for cell = 1:size(groupid,1)
    data = load(['F:\ClarkP30_Recordings\ProcessedData',filesep,groupid{cell,1}],'hdTuning');
    tuning = data.hdTuning{groupid{cell,5},groupid{cell,4}};
    smooth_tune = smooth_hd([tuning tuning tuning]);
    smooth_tune = smooth_tune(1, 61:120);
    mod_depth(cell,1) = (max(smooth_tune)+10e-20)/ (min(smooth_tune)+10e-20) ;
    
    tune(cell,:) = smooth_tune;
end

% Save data
shuffled_dic = table(groupid(:,1),groupid(:,2),groupid(:,3),groupid(:,4),groupid(:,5),z,p_val,shuff_pass,mod_depth,'VariableNames',{'session','tetrode','cell','baseline_index','cell_idx','z','p_val','shuff_pass','mod_depth'});
writetable(shuffled_dic,'d:/Users/BClarkLab/github/Berkowitz_et_al_2021/shuffled_dic_good.csv')


%% Plot distributions of all cells that meet critera
figure; 
edges = 0:.05:4;
histogram(hd_metrics.dic(shuff_pass == 1 & p_val < .001 & mod_depth > 1.5),edges)
hold on;
histogram(hd_metrics.dic(shuff_pass == 0),edges)


% %% 2. Examine shuffling results 
% df = readtable('d:/Users/BClarkLab/github/Berkowitz_et_al_2021/HD_Modulation/data/shuffled_mlv.csv');
% save_path = 'd:/Users/BClarkLab/github/Berkowitz_et_al_2021/Figures/HD_examples';
% 
% % keep cells that met for HD (v > .001 & z_score > 10 and passed 99th
% % percentile for shuffling). 
% df_keep = df(df.p_val < .01 & df.shuff_pass == 1 & df.z > 10,:);
% 
% %% Plot HD figs for each cell and save to folder
% for i = 1:length(df_keep.SessionID)
%     disp(['processing:', df_keep.SessionID{i}])
%     %Load data file and find < *Data*Session*Cell >
%     data = load(['/Users/lauraberkowitz/Documents/ProcessedData/',df_keep.SessionID{i}]);
%     
%     session = 1; % look at session 1
%     tet = sscanf(df_keep.tetrode{i},'TT%d.mat');
%     cell = find_cells(data,tet,df_keep.cell(i));
%     
%     HD_plots(data,session,cell,save_path,0)
% end

function smoothed = smooth_hd(tuning)
gaus = @(tuning,mu,sig)exp(-(((tuning-mu).^2)/(2*sig.^2)));
window = 10;
x = -window/2:window/2;
mu = 0;
sig = 3;
kernel = gaus(x,mu,sig);
smoothed = conv(tuning,kernel/sum(kernel),'same');
end
