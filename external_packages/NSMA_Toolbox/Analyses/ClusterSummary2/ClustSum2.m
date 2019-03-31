function [fh] = ClustSum2(fn, T, X, Y, epochs, iCell, epochs_plotted)

% ClustSum2  Produces sheet with: log ISI vs. time; plots of firing locations
%
% [fh] = ClustSum2(fn, T, X, Y, epochs, iCell)
%
% INPUTS:
%    fn - filename of cluster file
%    T - ts of timestamps for a single cluster%        
%    epochs - a struct object with fields:
%       epochs.names = cell array of strings of epoch names 
%       epochs.intervals = cell array of 1x2 arrays with [start_ts  end_ts] (start
%           and end timestamps of each epoch) -- elements in cell array correspond
%           to epochs, listed in same sequence as in epochs.names
%    X,Y - tsds of x,y position (pixel) data
%    iCell - cell number (in tfiles list for that session)
%    epochs_plotted (optional) - names of epochs (must match those in epochs struct
%       described above) to show in scatterplots of spike locations.  A cell array
%       of similar form to epochs.names
%
% OUTPUTS: 
%   fh - the figure handle, so that this function can be run from an overall program 
%       that loops over all the cells and produces multiple graphs
%
% PL 2002, last modified '04 by MN


% PREPROCESSING

% convert underscores to dashes
fn = strrep(fn,'_','-');
fn = strrep(fn,'/','//');

% find start and end times
EndingTime = EndTime(T);
StartingTime = StartTime(T);
   
% create figure with respect to screensize
sc = [0 0 1 1];
fh = figure(1);
set(fh,'Units','normalized','Position',[sc(1)+.05,sc(2)+.1,sc(3)*0.9,sc(4)*0.8]);
orient tall;
    
set(gcf, 'Name', fn, 'NumberTitle', 'Off');
clf;


% PLOT: LogISI vs time
subplot(3, 1, 1);
t = Data(T);

if ~isempty(t)
    logISI = log10(diff(t));
    plot(t(1:end-1)/(60*10000),logISI,'.','MarkerSize', .5);  
    y0 = max(logISI);
    y1 = min(logISI); 
else    % if no spikes, draw a plot anyway showing the entire timespan with no spikes in it
    t = [epochs.intervals{1}(1) epochs.intervals{end}(2)];
    logISI = [2 6];
    plot(t/(60*10000),logISI,'-w');
    y0 = 6; y1 = 2;
end

%axis tight;

%xtitle = mean([epochs.intervals{1}(1)/(60*10000)  epochs.intervals{1}(2)/(60*10000)]);
xtitle = epochs.intervals{1}(2)/(60*10000);
%text(xtitle,y0+(y0-y1)/5.5,fn);
text(.125,1+1/5.5,fn,'Units','normalized');
%title(fn);
ylabel('log10(ISI)');
xlabel('min');
axis([epochs.intervals{1}(1)/(60*10000)-.15 epochs.intervals{end}(2)/(60*10000)+.15 y1 y0]);
hold on;   
DrawEpochLines(y0,y1,epochs);
hold off;
drawnow


% PLOT: ScatterPlots of spike locations  

if nargin == 7
    str = epochs_plotted{1};
    start_ts = [];
    end_ts = [];
    j = 1;
    for i=1:length(epochs.names)
        if j<=length(epochs_plotted)
            if strcmp(epochs.names{i},epochs_plotted{j}) 
            %if epochs.names{i} == epochs_plotted{j} 
                if j>1
                    str = [str ', ' epochs_plotted{j}];
                end %if
                start_ts = [start_ts; epochs.intervals{i}(1)];
                end_ts = [end_ts; epochs.intervals{i}(2)];
                j = j+1;
            end %if
        end%if
    end %for
    T = restrict(T,start_ts,end_ts);
end %if
                                       
% Get location data for where spikes & stims occurred:
[XSpikes, YSpikes] = ScatterFields({T}, X, Y);
                                                            
% draw a plot with spikes for each cell:
axes('position',[0.2 0.2 0.3 0.3]);
%subplot(nPlot, 1, 2);
plot(Data(XSpikes), Data(YSpikes), '.');

axis equal;
axis tight;
axis([min(Data(X)) max(Data(X)) min(Data(Y)) max(Data(Y))]);
set(gca,'XTick',[]); set(gca,'XTickLabel',{});
set(gca,'YTick',[]); set(gca,'YTickLabel',{})
title('Locations of spikes');
if nargin == 7
    text(.6,-.08,['epochs shown in these 2 plots: ' str],'Units','normalized');
end %if
drawnow;


% PLOT: PlaceField (x by y) (Occupancy normalized!)

YtoXratio = (max(Data(Y))-min(Data(Y)))/(max(Data(X))-min(Data(X)));
axes('position',[0.5 0.2 0.3 0.3]);
%subplot(nPlot, 1, 3);
axis off
[PF, Occ] = TuningCurves(T,X,64,Y,64*YtoXratio);
f = find(PF > 0 & Occ == 0); PF(f) = 0;
f = find(PF == NaN); PF(f) = 0;
PF = PF ./ (Occ + 0.0001);
PF = flipud(PF);
[npf,mpf] = size(PF);
PFvals = reshape(PF,1,npf*mpf);

%zmax = prctile(reshape(PF,1,npf*mpf),99.5)
    
% write own percentile function because # users of Matlab Stats. Toolbox is limited &
% above command sometimes causes it to bomb
PFsort = sort(PFvals);
x = [floor(.995*length(PFsort)) floor(.995*length(PFsort))+1];
y = PFsort(x);
zmax = interp1(x,y,.995*length(PFsort));
    
imagesc(PF, [0 max(zmax,5)]);
colormap(1-hot);
        
axis equal;
axis tight;
set(gca, 'XTick', [], 'YTick', []);
title('Place Field (Occ-normd)')
hold off
drawnow



%--------------------------------------------
function DrawEpochLines(y0,y1,epochs)

    for iep = 1:length(epochs.names)
        t0 = epochs.intervals{iep}(1)/(60*10000);
        t1 = epochs.intervals{iep}(2)/(60*10000);
        tm = (t0+t1)/2;
        hhl = line([t0 t0], [y0 y1]);
        set(hhl,'Color','r');
        hhl = line([t1 t1], [y0 y1]);
        set(hhl,'Color','b');
        text(tm,y0-(y1-y0)/15,epochs.names{iep});
    end
