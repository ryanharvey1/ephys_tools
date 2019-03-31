function EEG1 = Filter4Theta(EEG0, StartTime, EndTime)
%
% Filter4Theta
%    EEGout = Filter4Theta(EEGin, StartTime, EndTime)
%
% INPUTS
%        EEGin = [c]tsd of eeg data
%        StartTime, EndTime = timerange within which to filter data
%
% OUTPUTS
%        EEGout = [c]tsd of eeg data filtered "appropriately"
%    
% ALGO
% Uses a YuleWalker linearly constructed filter to filter 
% the EEG signal between 6 and 10 Hz.
% 

% ADR 1998 
% version L4.0
% Status PROMOTED

dt = DT(EEG0);
sampling_frequency = 10000/dt;
Nyq = sampling_frequency/2; 		% Nyquist frequency

R = Restrict(EEG0, StartTime - (dt * 2^10), EndTime + (dt * 2^10));
%R = InterpolateEEG(R);

f = [0 4/Nyq 6/Nyq 10/Nyq 12/Nyq 1];
m = [0 0     1     1      0      0];
S = warning;
warning off
[b,a] = yulewalk(50, f, m);
warning(S);
%[h,w] = freqz(b, a);
%plot(f,m,w/pi,abs(h));


EEGd = filtfilt(b, a, Data(R));
switch class(EEG0)
case 'tsd'
   EEG1 = tsd(Range(R, 'ts'), EEGd);
case 'ctsd'
   EEG1 = ctsd(Start(R), DT(R), EEGd);
otherwise
   error('Unknown EEG input class.');
end

EEG1 = Restrict(EEG1, StartTime, EndTime);

