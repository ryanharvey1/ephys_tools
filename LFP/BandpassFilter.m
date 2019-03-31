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
% taken from CMBHOME
      

Wn_theta = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; % normalized by the nyquist frequency

[btheta,atheta] = butter(3,Wn_theta);

signal_filtered = filtfilt(btheta,atheta,signal);