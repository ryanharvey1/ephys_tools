%% Import data
dataset = 'F:\ClarkP30_Recordings\ProcessedData\';
save_path = 'F:\ClarkP30_Recordings\Analysis\';
hd = readtable('F:\ClarkP30_Recordings\Analysis\hd_cell_list.csv');
cell = readtable('F:\ClarkP30_Recordings\Analysis\hd_cell_spike_idx.csv');
cell = cell.n;
cell_class = readtable('F:\ClarkP30_Recordings\Analysis\Cell_Classification\processed\pyr_int_df.csv');
cd(dataset)

mean_lin_vel = [];
mean_ang_vel=[];
total_distance = [];

peak_rate = [];
mean_rate = [];
nSpikes = [];
files = unique(hd.SessionID);
%% Loop through data
for i = 1:length(files)
    % Load data structure
    data = load([dataset,files{i}],'measures','BasicLoco','varnames');
    % Average Linear Velocity
    mean_lin_vel(i,1) = data.BasicLoco.MeanVelocity(1,1);
    % Average Angular Velocity
    mean_ang_vel(i,1) = data.BasicLoco.AverageAnglePerSec(1,1);
    % Movement
    total_distance(i,1) = data.BasicLoco.OverallDistTraveled(1,1);
    
    % Peak Firing Rate (from measures)
    peak_rate = [peak_rate; data.measures(:,contains(data.varnames,'PeakRate'),1)];
    % Average Firing Rate (from measures)
    mean_rate = [mean_rate; data.measures(:,contains(data.varnames,'OverallFiringRate'),1)];
    % Total nSpikes (to remove those with <100
    nSpikes = [nSpikes; data.measures(:,contains(data.varnames,'nSpikes'),1)];
end


