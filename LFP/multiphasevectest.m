

% testpower_phase
% thetarange=4:0.25:12;
% gammarange=30:1:120;


function [phase,pow]=multiphasevectest(f,S,Fs,width)
% FUNCTION [phase,pow]=multiphasevec(f,S,Fs,width)
%
% Returns the phase and power as a function of time for a range of
% frequencies (f).
%
% Simply calls phasevec in a loop.
%
% INPUT ARGS:
%   f = [2 4 8];   % Frequencies of interest
%   S = dat;       % Signal to process
%   Fs = 256;      % Sampling frequency
%   width = 6;     % Width of Morlet wavelet (>= 5 suggested).
%
% OUTPUT ARGS:
%   phase- Phase data [freqs,time]
%   power- Power data [freqs,time]
%

pow = nan(length(f),length(S)); 
phase = nan(length(f),length(S));

for a=1:length(f)
     [phase(a,:),pow(a,:)]=phasevec(f(a),S,Fs,width);
end

end

%%
function thetaPhs = extractThetaPhase(root,method,band)

if isempty(root.lfp)
	error('Please load lfp with active_lfp before calling');
end

switch(lower(method))
	
	case 'waveform'
		if exist('band','var') || ~isempty(band)
			warning('band not used when computing phase by waveform')
		end
		
		%  method described by 
		broadBand = [1 60];
		filtSig = buttfilt(root.lfp.signal,broadBand,root.lfp.fs,'bandpass',4);
		
		% define window detection params
		min_sep = root.lfp.fs / 20; % 50ms separation between peak and next trough
		min_length = min_sep; % same for the separation between trough and next peak
		
		% Find troughs and peaks 
		% as onsets and offsets of signal with positive slope
		trphsAndPeaks = CMBHOME.Utils.OverThresholdDetect(sign(diff(filtSig)),0,min_sep,min_length);
		
		% Find zero crossings as sign changes in signal
		upAndDownZeroXings = CMBHOME.Utils.OverThresholdDetect(sign(filtSig),0,min_sep,min_length);

		waveLandmarks = sort([trphsAndPeaks, upAndDownZeroXings]);
		
		error('This code is not complete')
		% this method was to be based off of the methods described by Belluscio
		% et al 2012, however they had a large number of channels
		% simultaneously recorded from which they computed the median signal
		% and then computed the wavefrom from this.  This likely dramatically
		% reduces the number of local minima and maxima found and thereby makes
		% this a practical endeavor.  To make this work on a signal channel of
		% data, the signal would likely have to be filtered more thereby
		% reducing the value of using this method.
		
	case 'wavelet'
		% returns a vector of values, spanning theta range, per time point 
		if exist('band','var') & ~isempty(band)
			thetarange = band;
		else
			thetarange = 4:0.1:12;
		end
		[thetaPhs,~] = multiphasevec(thetarange,root.lfp.signal,root.lfp.fs,16);

		
	case 'hilbert'
		% returns a single phase estimate by time point using bandpass filtering
		if exist('band','var') & ~isempty(band)
			thetarange = band;
		else
			thetarange = [4 12];
		end
		thetaPhs = angle(hilbert(buttfilt(root.lfp.signal,thetarange,root.lfp.fs,'bandpass',4)));
	
	otherwise
		error('Unknown method')
end
end

%%
function y = morlet(f,t,width)
% function y = morlet(f,t,width)
% 
% Morlet's wavelet for frequency f and time t. 
% The wavelet will be normalized so the total energy is 1.
% width defines the ``width'' of the wavelet. 
% A value >= 5 is suggested.
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)
%
% See also: PHASEGRAM, PHASEVEC, WAVEGRAM, ENERGY 
%
% Ole Jensen, August 1998 

sf = f/width;
st = 1/(2*pi*sf);
A = 1/sqrt(st*sqrt(pi));
y = A*exp(-t.^2/(2*st^2)).*exp(i*2*pi*f.*t);
end

%%
function [y, pow] = phasevec(f,s,Fs,width)
% FUNCTION [y,pow] = phasevec(f,s,Fs,width)
%
% Return the phase as a function of time for frequency f. 
% The phase is calculated using Morlet's wavelets. 
%
% INPUT ARGS:
%   f = 8;     % Frequency of interest
%   s = dat;   % The eeg signal to evaluate
%   Fs = 256;  % Sample rate
%   width = 6; % width of Morlet wavelet (>= 5 suggested).
%
% OUTPUT ARGS:
%   y- Phase as a function of time
%   pow- Power as a function of time
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);

t=-3.5*st:dt:3.5*st;
m = morlet(f,t,width);

y = conv(s,m);

pow = abs(y).^2;
pow = pow(ceil(length(m)/2):length(pow)-floor(length(m)/2));

l = find(abs(y) == 0); 
y(l) = 1;

% normalizes phase estimates to length one
y = y./abs(y);
y(l) = 0;
   
y = angle( y(ceil(length(m)/2):length(y)-floor(length(m)/2)) );
end
