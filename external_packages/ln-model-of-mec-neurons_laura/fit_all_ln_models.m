%% Description
% The model: r = exp(W*theta), where r is the predicted # of spikes, W is a
% matrix of one-hot vectors describing variable (P, H, S, or T) values, and
% theta is the learned vector of parameters.

%% compute the position, head direction, speed, theta phase, egocentric bearing, and distance matrices

% initialize the number of bins that position, head direction, speed, and
% theta phase will be divided into
n_pos_bins = 20;
n_dir_bins = 18;
n_speed_bins = 10;
n_theta_bins = 18;
n_ego_bins = 18;

% compute position matrix
[posgrid, posVec] = pos_map([posx_c posy_c], n_pos_bins, boxSize);

% compute head direction matrix
[hdgrid,hdVec,direction] = hd_map(hd,n_dir_bins);

% compute speed matrix
[speedgrid,speedVec,speed] = speed_map(posx_c,posy_c,n_speed_bins);

% compute theta matrix
[thetagrid,thetaVec,phase] = theta_map(filt_eeg,post,eeg_sample_rate,n_theta_bins);

% compute EBC matrices (bearing and distance)
[egogrid, ~, ego_bearing, distances, ego_bins, dist_bins,dist_range] = ebc_map(pos,hd, n_ego_bins);

% remove times when the animal ran > 100 cm/s (these data points may contain artifacts)
too_fast = find(speed >= 100);
posgrid(too_fast,:) = []; hdgrid(too_fast,:) = []; 
speedgrid(too_fast,:) = []; thetagrid(too_fast,:) = [];
spiketrain(too_fast) = []; %distgrid(too_fast,:) = [];
egogrid(too_fast,:) = [];


%% Fit all 15 LN models

numModels = 30;
testFit = cell(numModels,1);
trainFit = cell(numModels,1);
param = cell(numModels,1);
A = cell(numModels,1);
modelType = cell(numModels,1);

% ALL VARIABLES
A{1} = [ posgrid hdgrid egogrid thetagrid speedgrid]; modelType{1} = [1 1 1 1 1];
% FOUR VARIABLES
A{2} = [ posgrid hdgrid egogrid thetagrid]; modelType{2} = [1 1 1 1 0];
A{3} = [posgrid hdgrid egogrid speedgrid]; modelType{3} = [1 1 1 0 1];
A{4} = [posgrid hdgrid thetagrid speedgrid]; modelType{4} = [1 1 0 1 1];
A{5} = [posgrid egogrid thetagrid speedgrid]; modelType{5} = [1 0 1 1 1];
A{5} = [hdgrid egogrid thetagrid speedgrid]; modelType{6} = [0 1 1 1 1];

% THREE VARIABLES
A{6} = [egogrid thetagrid speedgrid]; modelType{6} = [0 0 1 1 1];
A{7} = [hdgrid thetagrid speedgrid]; modelType{7} = [0 1 0 1 1];
A{8} = [hdgrid egogrid speedgrid]; modelType{8} = [0 1 1 0 1];
A{9} = [hdgrid egogrid thetagrid]; modelType{9} = [0 1 1 1 0];
A{10} = [ posgrid thetagrid speedgrid]; modelType{10} = [1 0 0 1 1];
A{11} = [ posgrid  egogrid speedgrid]; modelType{11} = [1 0 1 0 1];
A{12} = [ posgrid egogrid thetagrid]; modelType{12} = [1 0 1 1 0];
A{13} = [ posgrid hdgrid speedgrid]; modelType{13} = [1 1 0 0 1];
A{14} = [ posgrid hdgrid thetagrid]; modelType{14} = [1 1 0 1 0];
A{15} = [ posgrid hdgrid egogrid]; modelType{15} = [1 1 1 0 0];

% TWO VARIABLE
A{16} = [ posgrid hdgrid]; modelType{16} = [1 1 0 0 0];
A{17} = [ posgrid egogrid ]; modelType{17} = [1 0 1 0 0];
A{18} = [ posgrid thetagrid]; modelType{18} = [1 0 0 1 0];
A{19} = [ posgrid speedgrid]; modelType{19} = [1 0 0 0 1];
A{20} = [ hdgrid egogrid]; modelType{20} = [0 1 1 0 0];
A{21} = [ hdgrid  thetagrid]; modelType{21} = [0 1 0 1 0];
A{22} = [ hdgrid speedgrid]; modelType{22} = [0 1 0 0 1];
A{23} = [ egogrid thetagrid]; modelType{23} = [0 0 1 1 0];
A{24} = [ egogrid speedgrid]; modelType{24} = [0 0 1 0 1];
A{25} = [ thetagrid speedgrid]; modelType{25} = [0 0 0 1 1];

% ONE VARIABLE 
A{26} = posgrid; modelType{26} = [1 0 0 0 0];
A{27} = hdgrid; modelType{27} = [0 1 0 0 0];
A{28} = egogrid; modelType{28} = [0 0 1 0 0];
A{29} = thetagrid; modelType{29} = [0 0 0 1 0];
A{30} = speedgrid; modelType{30} = [0 0 0 0 1];

% compute a filter, which will be used to smooth the firing rate
filter = gaussmf(-4:4,[2 0]); filter = filter/sum(filter); 
dt = post(3)-post(2); fr = spiketrain/dt;
smooth_fr = conv(fr,filter,'same');

% compute the number of folds we would like to do
numFolds = 10;

parfor n = 1:numModels
    fprintf('\t- Fitting model %d of %d\n', n, numModels);
    [testFit{n},trainFit{n},param{n}] = fit_model(A{n},dt,spiketrain,filter,modelType{n},numFolds);
end
