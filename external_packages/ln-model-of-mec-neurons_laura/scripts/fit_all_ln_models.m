%% Description
% The model: r = exp(W*theta), where r is the predicted # of spikes, W is a
% matrix of one-hot vectors describing variable (P, H, E, or S) values, and
% theta is the learned vector of parameters.
function [glm_res,testFit,trainFit,param,smooth_fr] = fit_all_ln_models(pos,hd,spiketrain,model_params,boxSize)
glm_res = struct;
glm_res.model_type = {'Place-HD-Ego';'Place-HD';'Place-Ego'; 'HD-Ego';'Place','HD','Ego'};

%% compute the position, head direction, speed, theta phase, egocentric bearing, and distance matrices

% initialize the number of bins that position, head direction, speed, and
% egocentric bearing will be divided into
n_pos_bins = model_params.pos_bins;
n_dir_bins = model_params.hd_bins;
n_speed_bins = model_params.speed_bins;
n_ego_bins = model_params.ego_bins;

% compute position matrix
[posgrid, ~] = pos_map(pos, n_pos_bins, boxSize);

% compute head direction matrix
[hdgrid,~,~] = hd_map(hd,n_dir_bins);

% compute speed matrix
[~,~,speed] = speed_map(pos(:,1),pos(:,2),n_speed_bins);

% compute EBC matrices (bearing and distance)
[egogrid, ~, ~, ~, ~, ~,~] = ebc_map(pos,hd, n_ego_bins);

% remove times when the animal ran > 100 cm/s (these data points may contain artifacts)
too_fast = find(speed >= 100);
posgrid(too_fast,:) = []; hdgrid(too_fast,:) = []; 
spiketrain(too_fast) = []; 
egogrid(too_fast,:) = [];


%% Fit all 7 LN models
% Initialize saved variables
testFit = cell(model_params.numModels,1);
trainFit = cell(model_params.numModels,1);
param = cell(model_params.numModels,1);
A = cell(model_params.numModels,1);
modelType = cell(model_params.numModels,1);

% ALL VARIABLES
A{1} = [ posgrid hdgrid egogrid]; modelType{1} = [1 1 1];

% TWO VARIBALE 
A{2} = [ posgrid hdgrid]; modelType{2} = [1 1 0];
A{3} = [ posgrid egogrid]; modelType{3} = [1 0 1];
A{4} = [ hdgrid egogrid]; modelType{4} = [0 1 1];

% ONE VARIABLE 
A{5} = posgrid; modelType{5} = [1 0 0];
A{6} = hdgrid; modelType{6} = [0 1 0];
A{7} = egogrid; modelType{7} = [0 0 1];

% compute a filter, which will be used to smooth the firing rate
filter = gaussmf(-4:4,[2 0]); filter = filter/sum(filter); 
dt = post(3)-post(2); fr = spiketrain/dt;
smooth_fr = conv(fr,filter,'same');

parfor n = 1:model_params.numModels
    fprintf('\t- Fitting model %d of %d\n', n, model_params.numModels);
    [testFit{n},trainFit{n},param{n}] = fit_model(A{n},dt,spiketrain,filter,modelType{n},model_params.numFolds);
end

% save model fits 
glm_res.train = cell2mat(trainFit);
glm_res.test = cell2mat(trainFit);
glm_res.param = cell2mat(param);

end



% Code for 5 predictor models

% ALL VARIABLES
% A{1} = [ posgrid hdgrid egogrid thetagrid speedgrid]; modelType{1} = [1 1 1 1 1];
% % FOUR VARIABLES
% A{2} = [ posgrid hdgrid egogrid thetagrid]; modelType{2} = [1 1 1 1 0];
% A{3} = [posgrid hdgrid egogrid speedgrid]; modelType{3} = [1 1 1 0 1];
% A{4} = [posgrid hdgrid thetagrid speedgrid]; modelType{4} = [1 1 0 1 1];
% A{5} = [posgrid egogrid thetagrid speedgrid]; modelType{5} = [1 0 1 1 1];
% A{5} = [hdgrid egogrid thetagrid speedgrid]; modelType{6} = [0 1 1 1 1];
% 
% % THREE VARIABLES
% A{6} = [egogrid thetagrid speedgrid]; modelType{6} = [0 0 1 1 1];
% A{7} = [hdgrid thetagrid speedgrid]; modelType{7} = [0 1 0 1 1];
% A{8} = [hdgrid egogrid speedgrid]; modelType{8} = [0 1 1 0 1];
% A{9} = [hdgrid egogrid thetagrid]; modelType{9} = [0 1 1 1 0];
% A{10} = [ posgrid thetagrid speedgrid]; modelType{10} = [1 0 0 1 1];
% A{11} = [ posgrid  egogrid speedgrid]; modelType{11} = [1 0 1 0 1];
% A{12} = [ posgrid egogrid thetagrid]; modelType{12} = [1 0 1 1 0];
% A{13} = [ posgrid hdgrid speedgrid]; modelType{13} = [1 1 0 0 1];
% A{14} = [ posgrid hdgrid thetagrid]; modelType{14} = [1 1 0 1 0];
% A{15} = [ posgrid hdgrid egogrid]; modelType{15} = [1 1 1 0 0];
% 
% % TWO VARIABLE
% A{16} = [ posgrid hdgrid]; modelType{16} = [1 1 0 0 0];
% A{17} = [ posgrid egogrid ]; modelType{17} = [1 0 1 0 0];
% A{18} = [ posgrid thetagrid]; modelType{18} = [1 0 0 1 0];
% A{19} = [ posgrid speedgrid]; modelType{19} = [1 0 0 0 1];
% A{20} = [ hdgrid egogrid]; modelType{20} = [0 1 1 0 0];
% A{21} = [ hdgrid  thetagrid]; modelType{21} = [0 1 0 1 0];
% A{22} = [ hdgrid speedgrid]; modelType{22} = [0 1 0 0 1];
% A{23} = [ egogrid thetagrid]; modelType{23} = [0 0 1 1 0];
% A{24} = [ egogrid speedgrid]; modelType{24} = [0 0 1 0 1];
% A{25} = [ thetagrid speedgrid]; modelType{25} = [0 0 0 1 1];
% 
% % ONE VARIABLE 
% A{26} = posgrid; modelType{26} = [1 0 0 0 0];
% A{27} = hdgrid; modelType{27} = [0 1 0 0 0];
% A{28} = egogrid; modelType{28} = [0 0 1 0 0];
% A{29} = thetagrid; modelType{29} = [0 0 0 1 0];
% A{30} = speedgrid; modelType{30} = [0 0 0 0 1];