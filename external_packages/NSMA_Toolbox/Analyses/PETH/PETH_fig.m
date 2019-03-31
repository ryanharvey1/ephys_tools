function PETH_fig(trigger_event_ts, time_half_window_msec, nbins, title_str, eeg_flag)

% PETH_fig  Produces PETH figure and saves as Matlab fig file
%
% PETH_fig(trigger_event_ts, time_half_window_msec, nbins, title_str, eeg_flag) 
%
% INPUTS:  
%   trigger_event_ts - timestamps of when events (e.g. stims) occurred
%   time_half_window_msec -  +/- time on either side of event (e.g. stim) that PETH should be calculated for
%   nbins - number of histogram bins for spike PETH
%   title_str - string to be prepended to figure title (e.g., rat #, session ID, other info.)
%   eeg_flag - if 'eeg' entered for this input, avg. eeg is superimposed on figure
% OUTPUTS:
%   (none)
%
% PL & MN 2002, last modified '04 by MN


% get spike data
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





% Produce fig
Ncells = length(S);
spikes_ts = [];
for icell = 1:length(S) 
     spikes_ts = sort([spikes_ts; Data(S{icell})]);
end

dt_msec = 2*time_half_window_msec/nbins;

sc = [0 0 1 1];
fh = figure;
set(fh,'Units','normalized','Position',[sc(1)+.05,sc(2)+.1,sc(3)*0.9,sc(4)*0.8]);
orient tall;

[M, SM, B, R] = PETH(spikes_ts, trigger_event_ts, dt_msec, nbins);
if Ncells>0
    M = M/Ncells;  SM = SM/Ncells;      % convert spikes/sec to spikes/cell/sec
end
SMabove = M+SM;  SMbelow = M-SM;
errors = reshape([SMbelow'; SMbelow'; SMbelow'; SMabove'; SMabove'; SMabove'; NaN*ones(1,nbins)], 7*nbins, 1);
binwidth=mean(diff(B));
bins = reshape([B-binwidth/3; B+binwidth/3; B; B; B-binwidth/3; B+binwidth/3; B], 7*nbins, 1);

if max(M) == 0
    ymax = 1;
else
    if max(SM) > 0
        ymax = max(errors);
    else
        ymax = max(M);
    end
end
sc1 = [-time_half_window_msec time_half_window_msec 0 ymax];

bar(B,M); hold on;
plot(bins,errors,'g'); 

axis(sc1);
%scale1=AXIS;
xlabel('time (ms)');
ylabel('avg. spikes/cell/sec');


% if 'eeg' option specified, get eeg data & add to graph
if nargin == 5
    if eeg_flag == 'eeg'
        
        eegfile=findfiles('CSC*');
        [d, fname, ext] = fileparts(eegfile{1});
        msec_eeg = time_half_window_msec*1.5;
        [avgEEG, errEEG, EEGts, EEGsamples] = PETH_EEG(trigger_event_ts, [fname ext], msec_eeg, msec_eeg);
        clear errEEG EEGsamples;
        
        AX1 = GCA;
        set(AX1,'box','off','YColor','b');
        AX2 = axes('position',get(GCA,'position'));
        plot(EEGts/10,avgEEG,'r');
        sc2=[-time_half_window_msec time_half_window_msec min(avgEEG) max(avgEEG)];
        axis(sc2);
        %scale2=AXIS;
        set(AX2,'box','on','color','none','YAxisLocation','right','xtick',[],'ytick',[],'YColor','r');
        set(get(AX2,'Ylabel'),'String','avg. EEG');
        
    end %if eeg_flag
end %if nargin


title(['PETH - ' title_str]);


% save figure as Matlab figure file
fnstrg = ['PETH_' title_str];
saveas(fh,[fnstrg '.fig']);
    
close all;