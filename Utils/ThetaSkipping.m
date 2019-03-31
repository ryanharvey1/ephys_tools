function R = ThetaSkipping(cor, lag,t_bin)
% Calculates the ratio of the difference between 1st and 2nd peaks in autocorrelogram
%
% See PAGE 996 in Sachin S. Deshmukh, D. Yoganarasimha, Horatiu Voicu and James J. Knierim
% J Neurophysiol 104:994-1006, 2010
%
% Adapted from Hasselmo toolbox
% Ryan H


cor(lag==0) = [];

lags = lag(round(end/2):end);

% filter signal between 1 and 10 hz
Wn_theta = [1/(t_bin^-1/2) 10/(t_bin^-1/2)]; % normalized by the nyquist frequency

[btheta,atheta] = butter(3,Wn_theta);

acs = filtfilt(btheta,atheta,cor);

acs = acs(round(end/2):end)+mean(cor);

[amp, ind] = extrema(acs);

l = lags(ind);

peak1 = max(amp(l>=.1 & l<=.2)); % first peak

peak2 = amp(find(l>.2, 1, 'first')); % second peak

if ~isempty(peak1) && ~isempty(peak2)
    R = (peak2-peak1) / max([peak1 peak2]); % theta skipping index
else
    R = NaN;
    disp('No peaks in theta range of autocorrelogram');
end

end