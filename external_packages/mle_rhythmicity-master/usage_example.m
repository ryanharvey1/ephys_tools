%% USAGE EXAMPLE
% This script will simulate a speed tuned oscillating cell and analyze the
% modulation of rhythmicity by running speed using rhythmicity_covar. 
%
% Copywrite (c) 2015,2016 Trustees of Boston University
% All rights reserved.
%
% This file is part of mle_rhythmicity revision 2.0. The last committed
% version of the previous revision is the SHA starting with 93862ac...
%
% This code has been freely distributed by the authors under the BSD
% license (http://opensource.org/licenses/BSD2-Clause). If used or
% modified, we would appreciate if you cited our paper:
%
% Climer JR, DiTullio R, Newman EL, Hasselmo ME, Eden UT. (2014),
% Examination of rhythmicity of extracellularly recorded neurons in the
% entorhinal cortex. Hippocampus, 25:460-473. doi: 10.1002/hipo.22383.

clear all;clc;close all;
%% Simulated cell parameters
SESSION_DURATION = 5*60;% Seconds
AVERAGE_SPEED = 30;% cm/sec
MAX_ACCELERATION = 30;% cm/sec^2, smooths speed to match
N = 100;% Search resolution for smooth
F_SNR = 50;% signal-to-noise ratio for speed to frequency
SHARP = 2;% Power to raise cosine to to sharpen rhythmicity tuning

Fs = 50;% Sample frequency, frames/second

mean_freq = 8;% Hz, mean frequency of rhythmicity
freq_range = 5;% Hz, the amount the frequency of the rhythmicity varies over the running speed
mean_fr = 6;% Hz, mean firing rate of the neuron
%% Simulate cell
data = struct;
data.ts = 0:1/Fs:SESSION_DURATION;% Time stamps for behavior
data.speed = exprnd(AVERAGE_SPEED,size(data.ts));

sim_ts = 0:1/(Fs*10):SESSION_DURATION;% Discreet bins for simulation
[~,i] = min(arrayfun(@(x)abs(mean(abs(diff(smooth(data.speed,x))*Fs))-MAX_ACCELERATION),exp(linspace(0,log(numel(data.ts)),N))));% Find the smooth that keeps acceleration reasonable
data.speed = smooth(data.speed,exp(i*log(numel(data.ts))/N));

% Model for example: sinusoidal cell locked to a sinusoid (with noise)
data.lambda = (cos(...
    cumsum(...
    interp1(...
    data.ts,max(awgn((data.speed-mean(data.speed))*freq_range/range(data.speed)+mean_freq,F_SNR),0),sim_ts,'spline')*2*pi/(Fs*10)))+1).^SHARP;
data.lambda = (mean_fr*SESSION_DURATION)/sum(data.lambda)*data.lambda;

% Generate random spike times under model
data.spike_ts = poissrnd(data.lambda);
data.spike_ts = arrayfun(@(i)(find(data.spike_ts>=i)-1)/(Fs*10),1:max(data.spike_ts),'UniformOutput',false);
data.spike_ts = cat(2,data.spike_ts{:});
data.spike_ts = sort(data.spike_ts+unifrnd(-1/Fs/20,1/Fs/20,size(data.spike_ts)));

% Interpolate to make speed vector
data.spike_speed = interp1(data.ts,data.speed,data.spike_ts);

%% Run analysis
[ params, confidence, stats, everything ] = rhythmicity_covar(...
    ... REQUIRED INPUTS
    data.spike_ts(:) ... Spike times
    ,SESSION_DURATION ... Duration of the recording
    ,data.spike_speed(:) ... Speed as a vector - this defines what may be used for parameters to covary with
    ...
    ... DEFINE WHICH PARAMETERS COVARY WITH WHICH COVARIATES
    ,'f_covar',[0 1] ... Defines the rhythmicity frequency as able to covary with the speed (1) and be offset by y-intercept (0)
    ,'tau_covar',[0 1] ... Defines the rhythmicity falloff as able to covary with the speed (1) and be offset by y-intercept (0)
    ,'a_covar',[0 1] ... Defines the rhythmicity falloff as able to covary with the speed (1) and be offset by y-intercept (0)
    ...
    ... OTHER STATISTICS TO RUN
    ,'post_hocs',struct('f',1,'a',1) ... Run post-hoc tests for the speed-frequency and speed-amplitude tuning
    );