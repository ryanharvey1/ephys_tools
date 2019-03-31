function EEG_PowerVsTime(ripple_thresh)

% EEG_PowerVsTime  Plots theta, gamma & ripple power and velocity vs. time; saves as Matlab .fig file
%
% EEG_PowerVsTime(ripple_thresh)
%
% INPUTS:
%   ripple_thresh - EEG envelope threshold (determined with FindRippleIntervals.m)
%       above which ripple intervals are defined to occur
% OUTPUTS:
%   (none)
% 
% Note: this program is very memory-expensive.  Expect bombs!  Needs some fixing to 
% alleviate this problem.
%
% MN 8/02, last modified '04 by MN


% Get eeg data
eegfiles=findfiles('CSC*');
[d, fname, ext] = fileparts(eegfiles{1});
eeg_tsd = ReadCR_tsd([fname ext]);


% get position data
if exist('position.mat')
    load position X Y;
end %if exist('position.mat')


% get epochs
if exist('epochs.mat')
    load epochs;
end %if exist('epochs.mat')


%Produce fig
tstart = StartTime(eeg_tsd);
tend = EndTime(eeg_tsd);
Xts = Range(X,'ts');
t0 = Xts(1);
t1 = Xts(end);

startts = max(tstart,t0);
endts = min(tend,t1);
xlims = [startts endts];
%clear tstart tend Xts t0 t1 startts endts;

fh = figure;
set(fh,'Units','normalized','Position',[.05,.1,0.9,0.8]);
orient tall;

subplot(4,1,1);
SpeedVsTime(X,Y,epochs,xlims);
set(gca, 'XTickLabel', []);
%text(.4,1.3,title_prefix,'Units','normalized');
%clear X Y; pack;

subplot(4,1,2);
ThetaPowerVsTime(eeg_tsd,5,epochs,xlims);
set(gca, 'XTickLabel', []);

subplot(4,1,3);
GammaPowerVsTime(eeg_tsd,5,epochs,xlims);
set(gca, 'XTickLabel', []);
%pack;

subplot(4,1,4);
RippleCountsVsTime(ripple_thresh,eeg_tsd,5,epochs,xlims);
%set(gca, 'XTickLabel', []);
%clear ripple_thresh; pack;


% Save figure as Matlab figure file
saveas(fh,'EEG_PowerVsTime.fig');
    
close all;