function [epochs, P_r, ts] = GammaEpochs_v2(self,active_lfp_ind, window, window_inc)
% Determines epochs where theta/delta power ratio is above criteria 
%
%
% INPUT: 
%   self: data.lfp
%
% OUTPUT: 
%   epochs: indices for 
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

epochs = []; % initialize
P_r = [];
ts = [];

if iscell(self.signal), disp('Only one epoch allowed'); return; end

sig_delta = DeltaFilter(self.signal(active_lfp_ind,:), self.lfpsamplerate);

sig_theta = BandpassFilter(self.signal(active_lfp_ind,:), self.lfpsamplerate, [5, 12]); %changed to 5-12hz to limit redundancy between frequency band delta (2-4hz)

delta_amplitude = InstAmplitude(sig_delta);

theta_amplitude = InstAmplitude(sig_theta);

s_i = 1:windowinc:length(self.signal(active_lfp_ind,:))-windowsamples+1;

ts = self.ts(round(windowsamples/2:windowinc:length(self.signal(active_lfp_ind,:))-windowsamples/2));

if length(ts)~=length(s_i), disp('the length of ts doesnt match s_i in thetaepochs'); keyboard; end

P_r = nan(length(s_i)-1, 1);

for i = 1:length(s_i)
   
    P_r(i) = mean(theta_amplitude(round(s_i(i):s_i(i)+windowsamples-1)).^2) / mean(delta_amplitude(round(s_i(i):s_i(i)+windowsamples-1)).^2);
    
end

P_r = theta_amplitude ./ delta_amplitude;

ind = OverThresholdDetect(P_r, thresh, merge_if*self.lfpsamplerate, remove_if*self.lfpsamplerate);   % returns all indeces for epochs longer than .5 seconds that meet threshold
                                                                                                    % from +CMBHOME/+Utils
if ~isempty(ind)
        
    epochs = self.ts(ind); % set epochs to valid running epochs
        
end

P_r = P_r(:);
ts = ts(:);
end

%###################### LOCAL FUNCTIONS ##################################

    function signal_filtered = DeltaFilter(signal, Fs)
        
        %% Takes 'signal' and bandpasses it for delta frequencies (2 to 4 Hz)
        
        % Arguments
        
        % signal - arbitrary signal to be filtered
        % Fs - sample frequency (specific to your recording rig) ex. Mark's
        % Neuralynx system has a sampling period of .000528, or 1893.9393... Hz
        
        Wn_theta = [2/(Fs/2) 4/(Fs/2)]; % normalized by the nyquist frequency
        
        [btheta,atheta] = butter(3,Wn_theta);
        
        signal_filtered = filtfilt(btheta,atheta,signal);
        
    end

    function signal_amplitude = InstAmplitude(signal)
        % Returns instantaneous amplitude using the Hilbert transform for 'signal'
        %
        % signal_phase = LFP.InstAmplitude(signal);
        
        signal_amplitude = abs(hilbert(signal));
        
    end

    function signal_filtered = BandpassFilter(signal, Fs, Fpass)
        % Takes 'signal' and bandpasses it to Fpass frequencies
        %
        % Arguments
        %
        % signal - arbitrary signal to be filtered
        % Fs - sample frequency
        % Fpass - 1x2 vector of F_low and F_high indicating passband
        %
        % signal_filtered = BandpassFilter(signal, Fs, Fpass)
        
        if exist('Fpass', 'var')
            if diff(Fpass)<=0
                help CMBHOME.LFP.BandpassFilter
                error('Fpass must be in format [f_low, f_high]');
            end
        else
            help CMBHOME.LFP.BandpassFilter
            error('You must pass argument Fpass. See help above.');
        end
        
        Wn_theta = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; % normalized by the nyquist frequency
        
        [btheta,atheta] = butter(3,Wn_theta);
        
        signal_filtered = filtfilt(btheta,atheta,signal);
    end

    function [ind,tf] = OverThresholdDetect(A, thresh, min_sep, min_length)
        % Searches for continuous epochs for which A is at least thresh, merges those
        % separated by less than min_sep indices, and returns ind, an Nx2 array of
        % indices for continuous epochs over thresh
        % ind = ThresholdDetect(A, thresh, min_sep, min_length);
        % [ind, epochs] = ThresholdDetect(A, thresh, min_sep, min_length);
        
        A = A(:);
        
        tmp_ind = find(A<thresh); % all A under thresh to be thrown out
        
        if all(A<thresh)
            
            ind = []; % none are over thresh
            
        elseif ~isempty(tmp_ind);
            
            ind(:,1) = [1;tmp_ind+1];
            
            ind(:,2) = [tmp_ind-1; length(A)];
            
            ind(ind(:,2)-ind(:,1) < 0, :) = []; % remove all contiguous crossings
            
            if isempty(ind), ind = [tmp_ind(1), tmp_ind(end)]; end
            
            % merge those close enough as per max_sep
            
            ind = reshape(ind', 1, numel(ind)); % reshape into vector like: start, stop, start2, stop2, start3...
            
            tmpdiff = diff(ind);
            tmpdiff(1:2:end) = min_sep+5; %make all inter-epoch times greater than max_sep
            
            mergeinds = find(tmpdiff<min_sep);
            
            mergeinds = cat(2, mergeinds, mergeinds+1);
            
            ind(mergeinds) = [];    % remove all overlapping epochs
            
            ind = reshape(ind, 2, numel(ind)/2)';
            
            % remove epochs less than min_length
            
            ind(ind(:,2) - ind(:,1)<min_length-1,:) = [];
            
        else
            
            ind = [1 length(A)];
            
        end
        
        
        tf = zeros(length(A),1); % build tf logical vector
        
        for i = 1:size(ind,1)
            
            tf(ind(i,1):ind(i,2)) = 1;
            
        end
        
        tf = logical(tf);
    end