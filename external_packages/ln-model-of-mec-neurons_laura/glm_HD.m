
function glm_res = glm_HD(data,session,celln,varargin)
%%% Description of run_me

% This script is segmented into several parts. 
% 1. data is loaded. 
% 2. HD LN model is fit to the cell's spike train. 
% The model uses head direction information to predict a section of the
% spike train
% 3. Model fitting and model performance is computed through
%    10-fold cross-validation, and the minimization procedure is carried out
%    through fminunc. 
% 4. Signrank test to determine if observed log-likelihood values differ
% from zero. If not, selected model is NaN, otherwise model selects for HD.
% 
% Input: 
%  - data: processed data file from ephys_tools. 
%  - session: session number dictates which event files are used to index
%  spikes/timestamps.
%  - celln: cell number used as index for data.Spikes. 
% Output:
%  - glm_res: structure that contains train/test data sets, parameters
%  for each LN model, and the selected model.  
%  - smooth_fr: Firing rate smoothed with guassian filter for plotting
%  purposes. 
%

% Code as implemented in Hardcastle, Maheswaranthan, Ganguli, Giocomo,
% Neuron 2017
% V1: Kiah Hardcastle, March 16, 2017
% Adapted for ephys_tools RH 2019 
% Added EBC predictors LB 02/2020
% Simplified for HD & made into function LB 09/2020 

p = inputParser;
addOptional(p,'numFolds',10,@isnumeric)
addOptional(p,'n_circ_bins',18,@isnumeric)
addOptional(p,'n_speed_bins',10,@isnumeric)
parse(p,varargin{:})

% Initialize model hyperparameters
model_params = struct;
model_params.numFolds = p.Results.numFolds;
model_params.hd_bins = p.Results.n_circ_bins;
model_params.speed_bins = p.Results.n_speed_bins;

% boxSize = length (in cm) of one side of the square box
boxSize=data.maze_size_cm(session);

frames=data.frames(data.frames(:,1)>data.events(1,session) &...
    data.frames(:,1)<data.events(2,session),:); 

% convert xy to cm 
frames(:,2)=rescale(frames(:,2),0,boxSize);
frames(:,3)=rescale(frames(:,3),0,boxSize);

% post = vector of time (seconds) at every 33 ms time bin
post=frames(:,1);

% spiketrain = vector of the # of spikes in each 33 ms time bin
spiketrain=histcounts(data.Spikes{celln}(data.Spikes{celln}>data.events(1,session) &...
    data.Spikes{celln}<data.events(2,session)),...
    linspace(post(1),post(end),length(post)+1))';

% posx_c = x-position in middle of LEDs
posx_c = frames(:,2);

% posy_c = y-position in middle of LEDs
posy_c = frames(:,3);
pos = [posx_c posy_c];

% head direction
hd = deg2rad(frames(:,4));

%% fit the model
fprintf('(2/5) Fitting all linear-nonlinear (LN) models\n')

[glm_res,smooth_fr] = fit_all_ln_models(pos,hd,spiketrain,post,model_params);


%% find the simplest model that best describes the spike train
fprintf('(3/5) Performing forward model selection\n')
selected_model = select_best_model(glm_res.test,model_params);

glm_res.best_model = selected_model; 
glm_res.smooth_fr = smooth_fr;

end


%% Local Functions 

function [glm_res,smooth_fr] = fit_all_ln_models(pos,hd,spiketrain,post,model_params)
%% Description
% The model: r = exp(W*theta), where r is the predicted # of spikes, W is a
% matrix of one-hot vectors describing variable (HD) values, and
% theta is the learned vector of parameters.
glm_res = struct;
glm_res.model_type = {'HD'};

%% Compute the head direction and speed

% initialize the number of bins that head direction and speed will be divided into
n_dir_bins = model_params.hd_bins;
n_speed_bins = model_params.speed_bins;

% compute head direction matrix
[hdgrid,~,~] = hd_map(hd,n_dir_bins);

% compute speed matrix
[~,~,speed] = speed_map(pos(:,1),pos(:,2),n_speed_bins);

% remove times when the animal ran > 100 cm/s (these data points may contain artifacts)
too_fast = find(speed >= 100);
hdgrid(too_fast,:) = []; 
spiketrain(too_fast) = []; 

%% Fit HD LN models
% HD model
A = hdgrid; 

% compute a filter, which will be used to smooth the firing rate
filter = gaussmf(-4:4,[2 0]); filter = filter/sum(filter); 
dt = mean(diff(post)); fr = spiketrain/dt;
smooth_fr = conv(fr,filter,'same');
num_folds = model_params.numFolds;

[testFit,trainFit,param] = fit_model(A,dt,spiketrain,filter,num_folds);

% save model fits 
glm_res.train = trainFit;
glm_res.test = testFit;
glm_res.param = param;

end

function [testFit,trainFit,param_mean] = fit_model(A,dt,spiketrain,filter,numFolds)

%% Description
% This code will section the data into 10 different portions. Each portion
% is drawn from across the entire recording session. It will then
% fit the model to 9 sections, and test the model performance on the
% remaining section. This procedure will be repeated 10 times, with all
% possible unique testing sections. The fraction of variance explained, the
% mean-squared error, the log-likelihood increase, and the mean square
% error will be computed for each test data set. In addition, the learned
% parameters will be recorded for each section.


%% Initialize matrices and section the data for k-fold cross-validation

[~,numCol] = size(A);
sections = numFolds*5;

% divide the data up into 5*num_folds pieces
edges = round(linspace(1,numel(spiketrain)+1,sections+1));

% initialize matrices
testFit = nan(numFolds,6); % var ex, correlation, llh increase, mse, # of spikes, length of test data
trainFit = nan(numFolds,6); % var ex, correlation, llh increase, mse, # of spikes, length of train data
paramMat = nan(numFolds,numCol);

%% perform k-fold cross validation
for k = 1:numFolds
    fprintf('\t\t- Cross validation fold %d of %d\n', k, numFolds);
    
    % get test data from edges - each test data chunk comes from entire session
    test_ind  = [edges(k):edges(k+1)-1 edges(k+numFolds):edges(k+numFolds+1)-1 ...
        edges(k+2*numFolds):edges(k+2*numFolds+1)-1 edges(k+3*numFolds):edges(k+3*numFolds+1)-1 ...
        edges(k+4*numFolds):edges(k+4*numFolds+1)-1]   ;
    
    test_spikes = spiketrain(test_ind); %test spiking
    smooth_spikes_test = conv(test_spikes,filter,'same'); %returns vector same size as original
    smooth_fr_test = smooth_spikes_test./dt;
    test_A = A(test_ind,:);
    
    % training data
    train_ind = setdiff(1:numel(spiketrain),test_ind);
    train_spikes = spiketrain(train_ind);
    smooth_spikes_train = conv(train_spikes,filter,'same'); %returns vector same size as original
    smooth_fr_train = smooth_spikes_train./dt;
    train_A = A(train_ind,:);
    
    opts = optimset('Gradobj','on','Display','off'); % removed ,'Hessian','on'
    
    data{1} = train_A; data{2} = train_spikes;
    if k == 1
        init_param = 1e-3*randn(numCol, 1);
    else
        init_param = param;
    end
    [param] = fminunc(@(param) ln_poisson_model(param,data),init_param,opts);
    
    %%%%%%%%%%%%% TEST DATA %%%%%%%%%%%%%%%%%%%%%%%
    % compute the firing rate
    fr_hat_test = exp(test_A * param)/dt;
    smooth_fr_hat_test = conv(fr_hat_test,filter,'same'); %returns vector same size as original
    
    % compare between test fr and model fr
    sse = sum((smooth_fr_hat_test-smooth_fr_test).^2);
    sst = sum((smooth_fr_test-mean(smooth_fr_test)).^2);
    varExplain_test = 1-(sse/sst);
    
    % compute correlation
    correlation_test = corr(smooth_fr_test,smooth_fr_hat_test,'type','Pearson');
    
    % compute llh increase from "mean firing rate model" - NO SMOOTHING
    r = exp(test_A * param); n = test_spikes; meanFR_test = nanmean(test_spikes); 
    
    log_llh_test_model = nansum(r-n.*log(r)+log(factorial(n)))/sum(n); %note: log(gamma(n+1)) will be unstable if n is large (which it isn't here)
    log_llh_test_mean = nansum(meanFR_test-n.*log(meanFR_test)+log(factorial(n)))/sum(n);
    log_llh_test = (-log_llh_test_model + log_llh_test_mean);
    log_llh_test = log(2)*log_llh_test;
    
    % compute MSE
    mse_test = nanmean((smooth_fr_hat_test-smooth_fr_test).^2);
    
    % fill in all the relevant values for the test fit cases
    testFit(k,:) = [varExplain_test correlation_test log_llh_test mse_test sum(n) numel(test_ind) ];
    
    %%%%%%%%%%%%% TRAINING DATA %%%%%%%%%%%%%%%%%%%%%%%
    % compute the firing rate
    fr_hat_train = exp(train_A * param)/dt;
    smooth_fr_hat_train = conv(fr_hat_train,filter,'same'); %returns vector same size as original
    
    % compare between test fr and model fr
    sse = sum((smooth_fr_hat_train-smooth_fr_train).^2);
    sst = sum((smooth_fr_train-mean(smooth_fr_train)).^2);
    varExplain_train = 1-(sse/sst);
    
    % compute correlation
    correlation_train = corr(smooth_fr_train,smooth_fr_hat_train,'type','Pearson');
    
    % compute log-likelihood
    r_train = exp(train_A * param); n_train = train_spikes; meanFR_train = nanmean(train_spikes);   
    log_llh_train_model = nansum(r_train-n_train.*log(r_train)+log(factorial(n_train)))/sum(n_train);
    log_llh_train_mean = nansum(meanFR_train-n_train.*log(meanFR_train)+log(factorial(n_train)))/sum(n_train);
    log_llh_train = (-log_llh_train_model + log_llh_train_mean);
    log_llh_train = log(2)*log_llh_train;
    
    % compute MSE
    mse_train = nanmean((smooth_fr_hat_train-smooth_fr_train).^2);
    
    trainFit(k,:) = [varExplain_train correlation_train log_llh_train mse_train sum(n_train) numel(train_ind)];

    % save the parameters
    paramMat(k,:) = param;

end

param_mean = nanmean(paramMat);

end

function [f, df] = ln_poisson_model(param,data)

%updated by LB to limit predictor to HD

X = data{1}; % subset of A
Y = data{2}; % number of spikes

% compute the firing rate
u = X * param;
rate = exp(u);

% roughness regularizer weight - note: these are tuned using the sum of f,
% and thus have decreasing influence with increasing amounts of data
b_hd = 5e1;

% start computing the Hessian
% rX = bsxfun(@times,rate,X);
% hessian_glm = rX'*X;

%% find the H parameters and compute their roughness penalties
% initialize parameter-relevant variables
J_hd = 0; J_hd_g = []; %J_hd_h = [];


% compute the contribution for f, df, and the hessian
[J_hd,J_hd_g,~] = rough_penalty_1d_circ(param,b_hd); % J_hd_h


%% compute f, the gradient, and the hessian

f = sum(rate-Y.*u) + J_hd;
df = real(X' * (rate - Y) + J_hd_g);
% hessian = hessian_glm + blkdiag(J_hd_h);
end

%% smoothing functions called in the above script
function [J,J_g,J_h] = rough_penalty_1d_circ(param,beta)

numParam = numel(param);
D1 = spdiags(ones(numParam,1)*[-1 1],0:1,numParam-1,numParam);
DD1 = D1'*D1;

% to correct the smoothing across first and last bin
DD1(1,:) = circshift(DD1(2,:),[0 -1]);
DD1(end,:) = circshift(DD1(end-1,:),[0 1]);

J = beta*0.5*param'*DD1*param;
J_g = beta*DD1*param;
J_h = beta*DD1;
end

function selected_model = select_best_model(testFit,model_params)
testFit_mat = testFit;
LLH_values = reshape(testFit_mat(:,3),model_params.numFolds,1);

% re-set if selected model is not above baseline
pval_baseline = signrank(LLH_values,[],'tail','right');

if pval_baseline > 0.05
    selected_model = NaN;
else
    selected_model = 1;
end
end
