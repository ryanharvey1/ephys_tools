% SpikesVsTime_script  Creates graph of spike counts vs. time
%
% Produces and saves "Spikes vs. Time" graph for a cell or all cells in a
%   particular group combined.
% Run from folder containing "TT" subdirectories (.t files within these subdirs).
% Make sure 'epochs.mat' saved workspace is in this folder, or "epochs" 
%   variable is in current workspace. 
% This script should also be in this folder, or in a folder in your path.
% The function 'DrawEpochLines' should be in your path.
% Variables S (cell array of ts objects of spikes in each cluster) & tfiles (cell
%   array of t-file names) retained in workspace in case you want to save them.
%
% PL,MN 2002, last modified '04 by MN


% get epochs
if exist('epochs.mat')
    load epochs;
end %if exist('epochs.mat')


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


% Produce figure

binsize_msec = 200;  % <- alter this to whatever binsize you want

Ncells = length(S);
spikelist = [];
for icell = 1:Ncells 
     spikelist = sort([spikelist; Data(S{icell})]);
end

sc = [0 0 1 1];
fh = figure(1);
set(fh,'Units','normalized','Position',[sc(1)+.05,sc(2)+.1,sc(3)*0.9,sc(4)*0.8]);
orient tall;

firstbin = round(epochs.intervals{1}(1)/10/binsize_msec);
lastbin = round(epochs.intervals{end}(2)/10/binsize_msec);

bincenters = [firstbin:lastbin];
counts = hist(spikelist/10/binsize_msec, bincenters);  %spikes/bin
counts = counts/Ncells/(binsize_msec/1000);  % spikes/cell/sec
bincentersMinutes = bincenters*binsize_msec/1000/60;

xstart = bincentersMinutes(1);
xend = bincentersMinutes(end);
xquarter = xstart+(xend-xstart)/4;
xhalf = xstart+(xend-xstart)/2;
x3quarters = xend-(xend-xstart)/4;
[newmax, i] = max(counts);
if i==1
    sortcounts = sort(counts);
    newmax = sortcounts(end-1);
end %if

subplot(4,1,1);
H = plot(bincentersMinutes,counts);
set(H(1), 'Color', [0 0.9 0]);
y0 = newmax; y1 = 0;
hold on; DrawEpochLines(y0,y1,epochs);
axis([xstart-.5*binsize_msec/60000 xquarter 0 newmax]);
text(xstart+(xquarter-xstart)/3.5,y0-(y1-y0)/5,['Spike Counts vs. Time -- total of ' num2str(length(S)) ' cells']);

subplot(4,1,2);
H = plot(bincentersMinutes,counts);
set(H(1), 'Color', [0 0.9 0]);
y0 = newmax; y1 = 0;
hold on; DrawEpochLines(y0,y1,epochs);
axis([xquarter xhalf 0 newmax]);

subplot(4,1,3);
H1 = plot(bincentersMinutes,counts);
set(H1(1), 'Color', [0 0.9 0]);
y0 = newmax; y1 = 0;
hold on; DrawEpochLines(y0,y1,epochs);
axis([xhalf x3quarters 0 newmax]);
text(-0.07,0.5,'avg.spikes/cell/sec','Units','normalized','Rotation',90);

subplot(4,1,4);
H1 = plot(bincentersMinutes,counts);
set(H1(1), 'Color', [0 0.9 0]);
y0 = newmax; y1 = 0;
hold on; DrawEpochLines(y0,y1,epochs);
axis([x3quarters xend+.5*binsize_msec/60000 0 newmax]);
xlabel('time (min)');


% Save as Matlab fig file
saveas(fh,'SpikesVsTime.fig');

close all;

clear TTfolders i tfiles_temp S_temp binsize_msec Ncells spikelist icell sc fh firstbin lastbin ...
    bincenters* counts xstart xend xquarter xhalf x3quarters sortcounts newmax H y0 y1;


