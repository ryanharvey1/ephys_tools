function [phase,pow]=multiphasevec(f,S,Fs,width)
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