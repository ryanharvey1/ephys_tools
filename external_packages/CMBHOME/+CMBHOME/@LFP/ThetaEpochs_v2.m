function [epochs, P_r, ts] = ThetaEpochs_v2(self, window, window_inc)
% Determines epochs where theta/delta power ratio is above criteria 
%


% Returns continuous theta/delta power ratio using instantaneous power from
% Hilbert transform (P_r). moving window average of each instant. power. vector. 
% Note that the first windows length worth of data will have edge effects,
% thus only the data - windowsize (number of samples in window (sec)) is
% returned in P_r.
%
% Epochs for which (P_r) is greater than 2 are
% stored in 'epochs' [onset, offset;...]. Epochs closer than 50 ms are merged.
% Epochs less than 10 seconds are removed.
%
% [epochs, P_r, ts] = root.lfp.ThetaEpochs(windowsize);

% Modified to accomodate B.ClarkLab data structure LB 9/2018

thresh = 2;

merge_if = .05; %(seconds) epochs which are this close to one another are merged

remove_if = 2.5; % (seconds) epochs less than this are removed

if ~exist('window', 'var')
    
    window = 1; % second
    
end

if ~exist('window_inc', 'var')
    
    window_inc = 1; % second
    
end

windowsamples = self.lfpsamplerate*window; % number of samples in window 
windowinc = self.lfpsamplerate*window_inc; % number of samples in window_inc

%smooth_kernel_std = 1; % seconds

%kernel = pdf('norm', -5*smooth_kernel_std*self.fs:5*smooth_kernel_std*self.fs, 0, smooth_kernel_std*self.fs);

import CMBHOME.*

import CMBHOME.Utils.*

epochs = []; % initialize
P_r = [];
ts = [];

if iscell(self.signal), disp('Only one epoch allowed'); return; end

sig_delta = LFP.DeltaFilter(self.signal, self.lfpsamplerate);

sig_theta = LFP.BandpassFilter(self.signal, self.lfpsamplerate, [5, 12]); %changed to 5-12hz to limit redundancy between frequency band delta (2-4hz)

delta_amplitude = LFP.InstAmplitude(sig_delta);

theta_amplitude = LFP.InstAmplitude(sig_theta);

s_i = 1:windowinc:length(self.signal)-windowsamples+1;

ts = self.ts(round(windowsamples/2:windowinc:length(self.signal)-windowsamples/2));

if length(ts)~=length(s_i), disp('the length of ts doesnt match s_i in thetaepochs'); keyboard; end

P_r = nan(length(s_i)-1, 1);

for i = 1:length(s_i)
   
    P_r(i) = mean(theta_amplitude(round(s_i(i):s_i(i)+windowsamples-1)).^2) / mean(delta_amplitude(round(s_i(i):s_i(i)+windowsamples-1)).^2);
    
end

P_r = theta_amplitude ./ delta_amplitude;

ind = CMBHOME.Utils.OverThresholdDetect(P_r, thresh, merge_if*self.lfpsamplerate, remove_if*self.lfpsamplerate);   % returns all indeces for epochs longer than .5 seconds that meet threshold
                                                                                                    % from +CMBHOME/+Utils
if ~isempty(ind)
        
    epochs = self.ts(ind); % set epochs to valid running epochs
        
end

P_r = P_r(:);
ts = ts(:);