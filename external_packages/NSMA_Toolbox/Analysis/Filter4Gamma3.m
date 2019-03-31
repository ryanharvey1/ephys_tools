function [Y,Z] = Filter4Gamma2(IN_tsd, lowlimit, highlimit)
%
% Filter for Gamma EEG
% 
% INPUT:
%        IN_tsd     = The EEG data to be filtered (a tsd object).
%        lowlimit   = lower limit of the band pass filter 
%                     (default is 6 Hz)
%        highlimit  = upper limit of the band pass filter 
%                     (default is 10 Hz) Values above  Hz don't work
% 
% OUTPUT: 
%        Y        = tsd of filtered data that is subsampled.
%        Z        = tsd of filtered data that is NOT subsampled.
%
% function Y = Filter4Theta(IN_tsd, lowlimit, highlimit)

% cowen Wed Mar 24 14:08:41 1999
% MODIFIED: Maurer 9/14/2009
%WARNING: I REMOVED THE HIGH LIMIT ERROR IN ORDER TO FILTER FOR VALUES OVER
%100 Hz. USE AT YOUR OWN RISK.

sFreq = 400; % Hz. The desired sample frequency. It can't be too large for theta.
input_sFreq = 10000/mean(diff(Range(IN_tsd,'ts')))
interval = floor(input_sFreq/sFreq)
if interval == 0
   interval = 1;
end
crdata = Data(IN_tsd);
crtt = Range(IN_tsd,'ts');

if nargin == 1
  % Default band pass
  lowlimit  = 35;  % Hz
  highlimit = 80; % Hz
end
% if highlimit > 100
%   error('Keep the upper limit below 100. The filter cannot handle larger values')
% end


% Filter parameters
F_Ny = sFreq/2;           % Hz
lowcut = lowlimit/F_Ny;   % Hz
highcut = highlimit/F_Ny; % Hz
N = 4;             % Order of the filter
passband = [lowcut highcut];  
ripple = .1;

% Make sure the Cr tsd is OK.
if length(crdata) ~= length(crtt)
  error('input is misaligned');
end



%[Bb,Ab] = butter(N, passband);
[Bc,Ac] = cheby1(N, ripple, passband);
%h = [abs(hh) abs(freqz(Bb,Ab,n)) abs(freqz(Bc,Ac,n))];
% Use filtfilt instead of filter because it corrects for phase
% distortion due to filtering. It runs through the filter twice.
F = filtfilt(Bc,Ac,crdata(1:interval:end));

% The interval is used for subsampling the data
Y = tsd(crtt(1:interval:end),F);

% upsample the filtered signal to original sampling frequency 
% and adjust length to original CR-data length
z = interp(F,interval);
size(z)
lcr = length(crdata)
lz = length(z);
if lz < lcr
   z = [z; zeros(lcr-lz,1)];
end
if lz > lcr
   z = z(1:lcr);
end
size(z)
Z = tsd(crtt,z);
