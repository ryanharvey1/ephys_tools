function SpikesVsThetaPhase(ts_MazeStart, ts_MazeEnd)


% SpikesVsThetaPhase  Produces figure of spike counts over theta cycle, saves as Matlab fig file
%
% SpikesVsThetaPhase(ts_MazeStart, ts_MazeEnd) 
%
% INPUTS:  
%   ts_MazeStart, ts_MazeEnd - start & stop timestamps of desired behavior period.
% OUTPUTS:
%   (none)
%
% If no inputs specified, program searches for "epochs.mat" workspace containing "epochs" struct 
% variable with start & stop timestamps of all epochs in session.  Assumes that the 1st & last
% epochs in this struct are the pre- & post- behavior rest periods.  If not the case, function 
% can be edited.
%
% MN 2002, last modified '04 by MN


% get epochs if no start & stop timestamps input
if nargin < 1
    if exist('epochs.mat')
        load epochs;
    end %if exist('epochs.mat')
    ts_MazeStart = epochs.intervals{2}(1);     %or adjust if 1st & last epochs in "epochs" struct are
    ts_MazeEnd = epochs.intervals{end-1}(2);   % not the pre- & post-behavior rest periods
end %if nargin


% Get spike data
TTfolders = dir;
tfiles = {};
S = {};
for i=3:length(TTfolders)
    if TTfolders(i).isdir
        if ( strcmp(TTfolders(i).name(1:2),'TT') | strcmp(TTfolders(i).name(1:2),'tt') )

            cd(TTfolders(i).name)
            
            tfiles_temp = findfiles('*.t'); 
            S_temp = LoadSpikes(tfiles_temp);

            tfiles = [tfiles; tfiles_temp];
            S = [S; S_temp];
            
            cd ..
            
        end %if ( strcmp ...
    end %if TTfolders ...
end %for


% Get theta phase of spikes based on eeg data
eegfiles=findfiles('CSC*');
[d, fname, ext] = fileparts(eegfiles{1});
eeg_maze = ReadCR_tsd([fname ext], ts_MazeStart, ts_MazeEnd);
eeg_theta = Filter4Theta(eeg_maze,6,11);  clear eeg_maze;
PH = ThetaPhase(S,eeg_theta,ts_MazeStart,ts_MazeEnd);  clear eeg_theta;


% Produce figure
Ncells = length(PH);
phlist=[];
for i=1:Ncells
	phdata=data(PH{i});
    if ~isempty(phdata)
	    phlist=[phlist; phdata(:,1)];
    end
end

sc = [0 0 1 1];
fh = figure(1);
set(fh,'Units','normalized','Position',[sc(1)+.05,sc(2)+.1,sc(3)*0.9,sc(4)*0.8]);
orient tall;

bincenters = [.005:.01:.995];   % binsize = 100ths of a theta cycle. Alter if you want to change binsize.

counts = hist(phlist, bincenters);
counts = counts/sum(counts);  %normalize to proportion of counts in each bin
bincenters = [bincenters bincenters+1];
counts = [counts counts];
H = plot(bincenters,counts);
max=mean(counts)*2;

set(H(1), 'Color', [0 0.9 0]);
axis([0 2 0 max]);
ylabel('proportion of spikes -- binsize=.01 theta cycle');
title(['Histogram of spikes over theta cycle -- ' num2str(Ncells) ' cells']);
xlabel('theta phase');
set(gca,'XTickLabel',[[0:.2:1] [.2:.2:1]]);


% Save as Matlab fig file
saveas(fh,'SpikesVsThetaTime.fig');

close all;