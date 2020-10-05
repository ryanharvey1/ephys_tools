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
save_path = 'F:\ClarkP30_Recordings\Analysis\GLM\';
dataset = 'F:\ClarkP30_Recordings\ProcessedData\';

%% Assign rats to include in analysis (if not assigned, will default to all
% rats. 
 rats={'LB03','LB04','LB07','LB10'}; 
 transgenic = {'LB04','LB07','LB10'};
%% Run compile cell measures to obtain ID of all units with min 100 spikes &
% 1Hz peak firing rate. 
compile_cell_measures(rats,transgenic)


%% Runs through Hardcastle LN GLM for HD,Place, and HD/Place(Full model)
batch_hd_glm(dataset,save_path,'F:\ClarkP30_Recordings\Analysis\cell_list.csv')

%% For all cells classified as HD, find additional measures of interest,
% compile, and save as csv for analysis in R or Python. Saves copy to analysis folder and keeps data_all 
% (long format of all measures across max 4 conditions) in workspace. 
compile_HD_measures