
%% 1. Perform shuffling procedure on processed data 

% Load cell list 
groupid = readtable('/Users/lauraberkowitz/Documents/drift_analysis/cell_list.csv');
groupid = table2cell(groupid);

% Perform MVL shuffling on all cells
[shuff_pass,p_val,z]= shuff(groupid,'feature',{'mvl'},'nshuffle',1000);

% Save data
shuffled_mvl = table(groupid(:,1),groupid(:,2),groupid(:,3),z,p_val,shuff_pass,'VariableNames',{'SessionID','tetrode','cell','z','p_val','shuff_pass'});
writetable(shuffled_mvl,'/Users/lauraberkowitz/Documents/drift_analysis/shuffled_mlv.csv')

% Plot distributions of all cells that meet critera
figure; 
edges = -3:.5:25;
histogram(z(shuff_pass == 1 & p_val < .001),edges)
hold on;
histogram(z(shuff_pass == 0),edges)

% Plot WT vs Tg 
figure; 
histogram(z(shuff_pass == 1 & p_val < .001 & z > 3 & contains(groupid(:,1),'LB03')),edges)
hold on;
histogram(z(shuff_pass == 1 & p_val < .001 & z > 3 & ~contains(groupid(:,1),'LB03')),edges)

%% 2. Examine shuffling results 
df = readtable('/Users/lauraberkowitz/Documents/drift_analysis/shuffled_mlv.csv');
save_path = '/Users/lauraberkowitz/Documents/cell_classification/Figures/HD_Figs';

% keep cells that met for HD (v > .001 & z_score > 10 and passed 99th
% percentile for shuffling). 
df_keep = df(df.p_val < .01 & df.shuff_pass == 1 & df.z > 10,:);

%% Plot HD figs for each cell and save to folder
for i = 1:length(df_keep.SessionID)
    disp(['processing:', df_keep.SessionID{i}])
    %Load data file and find < *Data*Session*Cell >
    data = load(['/Users/lauraberkowitz/Documents/ProcessedData/',df_keep.SessionID{i}]);
    
    session = 1; % look at session 1
    tet = sscanf(df_keep.tetrode{i},'TT%d.mat');
    cell = find_cells(data,tet,df_keep.cell(i));
    
    HD_plots(data,session,cell,save_path,0)
end


