% Alzheimer_pipeline runs analysis scripts for TgF344-AD & Spatial Nav
% analysis LB August 2020
%   1. Initialize Paths (save path and processed data location)
%   2. Define rats to include (must be indicated by animal id in metadata)
%   3. Find cells with min 100 spikes and 1Hz peak firing rate
%   (compile_cell_measures)
%   4. Run LN GLM for head direction and shuffling procedures for
%   directional info content (batch_hd_glm)
%   5. Run extra measures for HD cell firing characteristics
%   (compile_HD_measures)

%% Add paths required for script
addpath('F:\ClarkP30_Recordings\AnimalMetaData');
addpath('F:\ClarkP30_Recordings\Analysis\')

%% Assign location to save data and to find data
dataset = 'F:\ClarkP30_Recordings\ProcessedData\';

%% Assign rats to include in analysis (if not assigned, will default to all
% rats. 
<<<<<<< HEAD
 rats={'LB03','LB04','LB07','LB10','PoS4'}; 
 transgenic = {'LB04','LB07','LB10'};
=======
 rats = {'LB03','LB04','LB07','LB10','LB11','PoS4','PoS3'}; 
 transgenic = {'LB04','LB07','LB10','PoS3'};
 
>>>>>>> b327fe90583f1fb0c83d8b12694bd8600adc68da
%% Run compile cell measures to obtain ID of all units with min 100 spikes &
% 1Hz peak firing rate. 
compile_cell_measures(rats,transgenic)

<<<<<<< HEAD
%% Import hd_cell_list and create index .csv using find cells (For python processing)
cell_list = readtable('/Users/lauraberkowitz/github/Berkowitz_et_al_2021/Cell_Classification/data/cell_list.csv');
cell_idx = [];
for i = 1 : length(cell_list.session)
    data = load(['/Users/lauraberkowitz/Google Drive/ClarkP30_Recordings/ProcessedData/',cell_list.session{i}],'spikesID');
    tet = sscanf(cell_list.tetrode{i},'TT%d.mat');
    cell_idx(i,1) = find_cells(data,tet,cell_list.cell(i));
=======
%% Import cell_list and create index .csv using find cells (For python processing)
cell_list = readtable('D:/Users/BClarkLab/github/Berkowitz_et_al_2021/Cell_Classification/data/cell_list.csv');
cell_list.cell_idx{1} = [];

for i = 1 : length(cell_list.session)
    data = load(['d:\Users\BClarkLab\GoogleDrive_\ClarkP30_Recordings\ProcessedData\',cell_list.session{i}],'spikesID');
    tet = sscanf(cell_list.tetrode{i},'TT%d.mat');
    cell_list.cell_idx{i,1} = find_cells(data,tet,cell_list.cell(i));
>>>>>>> b327fe90583f1fb0c83d8b12694bd8600adc68da
end
% convert cell to double
cell_list.cell_idx = cell2mat(cell_list.cell_idx);
% write back to list
writetable(cell_list,'D:/Users/BClarkLab/github/Berkowitz_et_al_2021/Cell_Classification/data/cell_list.csv')

%% Calculate metrics for cell classification using cell list
classify_pyr_int

%% Runs through Hardcastle LN GLM for HD,Place, and HD/Place(Full model)
% batch_hdVel_glm(dataset,save_path,'D:\Analysis\cell_list.csv')


<<<<<<< HEAD
cell_list.cell_idx = cell_idx;
writetable(cell_list,'/Users/lauraberkowitz/github/Berkowitz_et_al_2021/Cell_Classification/data/cell_list.csv')

%% Runs through Hardcastle LN GLM for HD,Place, and HD/Place(Full model)
% batch_hd_glm(dataset,save_path,'F:\ClarkP30_Recordings\Analysis\cell_list.csv')

=======
>>>>>>> b327fe90583f1fb0c83d8b12694bd8600adc68da
%% For all cells classified as HD, find additional measures of interest,
% compile, and save as csv for analysis in R or Python. Saves copy to analysis folder and keeps data_all 
% (long format of all measures across max 4 conditions) in workspace. 
compile_HD_measures