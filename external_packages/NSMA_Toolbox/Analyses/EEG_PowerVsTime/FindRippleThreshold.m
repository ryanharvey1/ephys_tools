function FindRippleThreshold(threshold_guess)

% FindRippleThreshold  Produces figures to visually check EEG envelope threshold for ripples
%
% FindRippleThreshold(threshold_guess)
%
% INPUTS:  
%   threshold_guess - your guess of EEG envelope threshold for ripples
% OUTPUTS:
%   (none)
%
% Produces a figure with 2 subplots:
%   1) Figure of the power envelope of ripple-filtered EEG vs. time, with visible threshold 
%       and intervals.  Solid black horizontal lines are plotted at the y-axis height of 
%       threshold_guess and span the x-axis range for each discrete ripple interval.
%   2) Ripple counts vs. time.  If threshold_guess is accurate, one would expect to see counts
%       of about 4-5 ripples per 5 sec. time bin.
%
% MN 4/04 (from PL code '00 & '01), last modified '04 by MN


% Get epochs
load epochs;


% Get eeg data
eegfiles=findfiles('CSC*');
[d, fname, ext] = fileparts(eegfiles{1});
eeg_tsd = ReadCR_tsd([fname ext]);


t0 = epochs.intervals{1}(1);
t1 = epochs.intervals{end}(end);
tstart = StartTime(eeg_tsd);
tend = EndTime(eeg_tsd);
startts = max(tstart,t0);
endts = min(tend,t1);
xlims = [startts endts];


% Filter for ripples
Z_ripple = Filter4Ripples(eeg_tsd,100,300);


% PRODUCE PLOTS

% Power env. vs. time with threshold & intervals shown
subplot(2,1,1);
[ripple_start_ts, ripple_end_ts] =  RippleIntervals(Z_ripple, 120, 1, epochs, xlims);

% Ripple counts vs. time
subplot(2,1,2);
RippleCounts(Z_ripple, ripple_start_ts, ripple_end_ts, epochs, xlims);
