%% Description of run_me

% This script is segmented into several parts. First, the data (an
% example cell) is loaded. Then, 15 LN models are fit to the
% cell's spike train. Each model uses information about 
% position, head direction, running speed, theta phase,
% or some combination thereof, to predict a section of the
% spike train. Model fitting and model performance is computed through
% 10-fold cross-validation, and the minimization procedure is carried out
% through fminunc. Next, a forward-search procedure is
% implemented to find the simplest 'best' model describing this spike
% train. Following this, the firing rate tuning curves are computed, and
% these - along with the model-derived response profiles and the model
% performance and results of the selection procedure are plotted.

% Code as implemented in Hardcastle, Maheswaranthan, Ganguli, Giocomo,
% Neuron 2017
% V1: Kiah Hardcastle, March 16, 2017
% Adapted for ephys_tools RH 2019 
% Added EBC predictors LB 02/2020

%% Clear the workspace and load the data

clear; clc

% load the data
fprintf('(1/5) Loading data from example cell \n')

data = load('LEM3206_S20190702114353.mat');
celln=find_cells(data,1,10);

postprocessFigures.main(data,...
    'cellid',{data.spikesID.TetrodeNum{celln},data.spikesID.CellNum(celln)},...
    'colorcode','HD');


session=3; %contains(data.mazetypes,'cylinder');

% description of variables included:
% boxSize = length (in cm) of one side of the square box
boxSize=data.maze_size_cm(session);

frames=data.frames(data.frames(:,1)>data.events(1,session) &...
    data.frames(:,1)<data.events(2,session),:); 
% 
% bad_frames=any(isnan(frames),2);
% 
% frames(bad_frames,:)=[];

% convert xy to cm 
frames(:,2)=rescale(frames(:,2),0,boxSize);
frames(:,3)=rescale(frames(:,3),0,boxSize);

% post = vector of time (seconds) at every 20 ms time bin
post=frames(:,1);

% spiketrain = vector of the # of spikes in each 20 ms time bin
spiketrain=histcounts(data.Spikes{celln}(data.Spikes{celln}>data.events(1,session) &...
    data.Spikes{celln}<data.events(2,session)),...
    linspace(post(1),post(end),length(post)+1))';

% posx_c = x-position in middle of LEDs
posx_c = frames(:,2);

% posy_c = y-position in middle of LEDs
posy_c = frames(:,3);

pos = [posx_c posy_c];

% filt_eeg = local field potential, filtered for theta frequency (4-12 Hz)
filt_eeg=data.lfp.theta(1,:)';

% eeg_sample_rate = sample rate of filt_eeg (250 Hz)
eeg_sample_rate=data.lfp.lfpsamplerate;

% sampleRate = sampling rate of neural data and behavioral variable (50Hz)
sampleRate=data.samplerate;

% head direction
hd=deg2rad(frames(:,4));

clear data frames

%% fit the model
fprintf('(2/5) Fitting all linear-nonlinear (LN) models\n')
tic
fit_all_ln_models
toc
%% find the simplest model that best describes the spike train
fprintf('(3/5) Performing forward model selection\n')
select_best_model

%% Compute the firing-rate tuning curves
fprintf('(4/5) Computing tuning curves\n')
compute_all_tuning_curves

%% plot the results
fprintf('(5/5) Plotting performance and parameters\n') 
plot_performance_and_parameters
