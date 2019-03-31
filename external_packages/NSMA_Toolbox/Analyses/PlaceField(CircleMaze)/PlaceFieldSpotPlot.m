function fh = PlaceFieldSpotPlot(V,cell_ts,horseshoe_flag,event_ts1,event_ts2,pf_center,titlestr,thetapow_tsd,varargin)

% PlaceFieldSpotPlot  Produces lap-by-lap raster of firing locations
% 
% fh = PlaceFieldSpotPlot(V,cell_ts,horseshoe_flag,event_ts1,event_ts2,pf_center,titlestr,thetapow_tsd,varargin)
% 
% INPUTS:
%   V - position tsd; angle on circular track in degrees [0..360) 
%   cell_ts - array of timestamps of spikes from one cell
%   horseshoe_flag - 1=animal runs in "horseshoe" pattern on circular track, where it runs once
%       around circle in one direction, turns and runs once around in reverse direction, etc.
%       0(or no inputs)=animal runs regular laps around circle
%   event_ts1 & 2 - array of timestamps of events for up to 2 different types of events.  These
%       inputs optional (may be empty [])
%   pf_center - location (degrees) of center of place field
%   titlestr - tetrode & cell no. to be used in fig title
%   thetapow_tsd - tsd of power at theta peaks, obtained from theta-filtered eeg
%   VARARGIN parameters:
%       1) may enter an array of timestamps, so that SpotPlot only consists of laps containing
%               those timestamps
%       2) may enter a "1" indicating that SpotPlot should consist of all laps NOT containing those
%               timestamps in first varargin
% OUTPUTS:
%   fh - figure handle
%
% Lap-by-lap raster showing positions of spikes and possibly other events.
%   Larger diameter of circles = higher firing rate.
%
% PL '01, last modified '04 by MN


% Count Laps  
% Check for time gaps in position data - if > 20 sec, make a break there and when data resumes
% after gap, count as a new lap
Vts = Range(V,'ts');
ixbeginepochs = [1; find(diff(Vts) > 20*10000)+1];
tsbeginepochs = Vts(ixbeginepochs);
ixendepochs = [find(diff(Vts) > 20*10000); length(Vts)];
tsendepochs = Vts(ixendepochs);
VtsTemp = [];  Vdata = [];  flagixlaps = [];  lastlapinepochs = [];  cumlaps = 0;
for i = 1:length(ixbeginepochs)
    if i<length(ixbeginepochs)
        Vepoch = restrict(V,Vts(ixbeginepochs(i)),Vts(ixbeginepochs(i+1)-1));
    else
        Vepoch = restrict(V,Vts(ixbeginepochs(i)),Vts(end));
    end %if
    [ix_laps_epoch, nLaps_epoch, V_epoch] = LapCount(Vepoch, horseshoe_flag);
    cumlaps = cumlaps + nLaps_epoch;
    lastlapinepochs = [lastlapinepochs; cumlaps];
    flagixlaps_epoch = zeros(size(data(Vepoch)));  
    flagixlaps_epoch(ix_laps_epoch) = 1;
    VtsTemp = [VtsTemp; Range(Vepoch,'ts')];  Vdata = [Vdata; data(Vepoch)];
    flagixlaps = [flagixlaps; flagixlaps_epoch];
end %for
V = tsd(VtsTemp,Vdata);
Vts = VtsTemp;  
ix_laps = find(flagixlaps == 1);
nLaps = length(ix_laps);
lapnums = 1:nLaps;

% Restrict all inputs to running periods
if 1  %if want to restrict everything to ALL running periods
    ts_end = tsendepochs;
end %if
if 0  % if want to restrict to only full (no partial) laps
    ts_end = Vts(ix_laps(lastlapinepochs));
end %if 
thetapow_tsd = restrict(thetapow_tsd,tsbeginepochs,ts_end);
cell_ts = data(restrict(ts(cell_ts),tsbeginepochs,ts_end));
event_ts1 = data(restrict(ts(event_ts1),tsbeginepochs,ts_end));
event_ts2 = data(restrict(ts(event_ts2),tsbeginepochs,ts_end));
flagixlaps = data(restrict(tsd(Vts,flagixlaps),tsbeginepochs,ts_end));
ix_laps = find(flagixlaps == 1);
V = restrict(V,tsbeginepochs,ts_end);
Vts = Range(V,'ts');  Vdata = data(V);


% Chop out laps that DO NOT (if only one varargin) or DO (if 2 varargins) contain 
%  any timestamps listed in 1st varargin
if length(varargin) > 0
    tslist = varargin{1};
    
    ixLast = 0;
    Vdelete = zeros(size(Vts));
    event1delete = zeros(size(event_ts1));
    event2delete = zeros(size(event_ts2));
    celldelete = zeros(size(cell_ts));
    thetapow_data = data(thetapow_tsd); thetapow_ts = Range(thetapow_tsd,'ts');
    thetapowdelete = zeros(size(thetapow_data));
    lapnums = 1:nLaps;
    flag_lapnums = zeros(size(lapnums));
    
    for i=1:nLaps
        ix1=ixLast;
        ix2=ix_laps(i);
        
        if length(varargin) == 1
            x = isempty(find(tslist >= Vts(ix1+1) & tslist <= Vts(ix2)));
        end % if 
        if length(varargin) == 2
            x = ~isempty(find(tslist >= Vts(ix1+1) & tslist <= Vts(ix2)));
        end %if 

        if x  %if lap is to be cut
            flag_lapnums(i) = 1;
            Vdelete(ix1+1:ix2) = 1;
            event1delete(find(event_ts1 >= Vts(ix1+1) & event_ts1 <= Vts(ix2))) = 1;
            event2delete(find(event_ts2 >= Vts(ix1+1) & event_ts2 <= Vts(ix2))) = 1;
            celldelete(find(cell_ts >= Vts(ix1+1) & cell_ts <= Vts(ix2))) = 1;
            thetapowdelete(find(thetapow_ts >= Vts(ix1+1) & thetapow_ts <= Vts(ix2))) = 1;
        end %if x
        
        ixLast = ix2;
    end %for i
    
    flag_laps = zeros(size(Vts));  flag_laps(ix_laps) = 1;
    
    % chop data flagged for deletion
    lapnums(find(flag_lapnums == 1)) = [];
    flag_laps(find(Vdelete == 1)) = [];
    Vts(find(Vdelete == 1)) = [];
    Vdata(find(Vdelete == 1)) = [];
    V = tsd(Vts,Vdata);
    event_ts1(find(event1delete == 1)) = [];
    event_ts2(find(event2delete == 1)) = [];
    cell_ts(find(celldelete == 1)) = [];
    thetapow_ts(find(thetapowdelete == 1)) = [];
    thetapow_data(find(thetapowdelete == 1)) = [];
    thetapow_tsd = tsd(thetapow_ts,thetapow_data);
    
    ix_laps = find(flag_laps == 1);
    nLaps = length(ix_laps);
        
end %if length(varargin) > 0;


%(previous version used INTERP1 function to get location data at spike/event times, but this screwed up at
% 0 -> 360 degree crossings)
V_at_cellts = data(restrict(V,cell_ts));
V_at_eventts1 = data(restrict(V,event_ts1));
V_at_eventts2 = data(restrict(V,event_ts2));

% caluclate Occupancy and Spike histos per Lap
LapSec = zeros(nLaps,1);
ixLast = 0;
for iL = 1:nLaps
    ts1 = Vts(ixLast+1);
    ts2 = Vts(ix_laps(iL));
    if horseshoe_flag == 1
        occcenters = -356:8:356;
    else
        occcenters = 2:4:358;
    end %if horseshoe
    [OccPerLap] = hist(Vdata(ixLast+1:ix_laps(iL)),occcenters);
    SV{iL} = V_at_cellts(  find(cell_ts >= ts1 & cell_ts <= ts2)  );
    eventV1{iL} = V_at_eventts1(  find(event_ts1 >= ts1 & event_ts1 <= ts2)  );
    eventV2{iL} = V_at_eventts2(  find(event_ts2 >= ts1 & event_ts2 <= ts2)  );
    [nnPerLap] = hist(SV{iL},occcenters);
    SpotSize{iL} = nnPerLap./(OccPerLap+1/10000);
    ix_zeros = find(SpotSize{iL} >= 9999);
    SpotSize{iL}(ix_zeros) = 0;        % set all bins to 0 where Occ is 0
    LapSec(iL) = (ts2-ts1)/10000;
    Vdata_lap = Vdata(find(Vts >= ts1 & Vts <= ts2));
    minv = min(Vdata_lap);
    maxv = max(Vdata_lap);
    if horseshoe_flag == 0
        Vnorun_ends{iL} = [0 minv(1) NaN maxv(1) 360];
    else
        Vnorun_ends{iL} = [-360 minv(1) NaN maxv(1) 360];
        if ~isempty(find(Vdata_lap>0))
            minv_mid = max(Vdata_lap(find(Vdata_lap<0)));
            maxv_mid = min(Vdata_lap(find(Vdata_lap>0)));
            Vnorun_mid{iL} = [minv_mid(1) maxv_mid(1)];
        else
            if iL == nLaps
                Vnorun_mid{iL} = [];
            end %if
        end %if ~isempty
    end %if horseshoe
    ixLast = ix_laps(iL);
end % for iL


sc = [0 0 1 1];
fh = figure;
set(fh,'Units','normalized','Position',[sc(1)+.05,sc(2)+.1,sc(3)*0.9,sc(4)*0.8]);
orient tall;


% PLOT LAP-SPIKE HISTOGRAM
axes('position', [0.1 0.87 0.8 0.08]);
lapSpikeCount = zeros(nLaps,1);
for ii=1:nLaps
    lapSpikeCount(ii) = length(SV{ii});    
end % for ii
bar(lapSpikeCount,'k'); hold on;

L = length(lapSpikeCount);
M = max(lapSpikeCount);

if M ~= 0;
    SpikesPerSec = lapSpikeCount./LapSec;
    SpikesPerSec_rescaled = SpikesPerSec*max(lapSpikeCount)/max(SpikesPerSec);
    plot(1:nLaps,SpikesPerSec_rescaled,'g'); %avg firing rate per lap
    LapSec_rescaled = LapSec*max(lapSpikeCount)/max(LapSec);
else
    LapSec_rescaled = LapSec/max(LapSec);
    M = 1;
    text(.42,.6,'(No Spikes)','Units','normalized','Color','r','FontSize',15);
end %if max...

plot(1:nLaps,LapSec_rescaled,'m'); %lap time in sec
  
if length(varargin) > 0
    set(gca, 'XTick', 1:nLaps, 'XTickLabel', lapnums);
end%if
axis([0  L+1  0  1.15*M]);
xlabel('Lap');
ylabel('Spikes');
title(titlestr);


% PLOT THETA POWER & VELOCITY VS. POSITION
axes('position', [0.1 0.72 0.8 0.08]);
if horseshoe_flag == 1
    ThetaPowerVsPosition(V, thetapow_tsd, -358:4:358, horseshoe_flag);  
else
    ThetaPowerVsPosition(V, thetapow_tsd, 1:2:359, horseshoe_flag);  
end %if horseshoe
hold on;
ylim = get(gca,'YLim');  ymax = ylim(2);
line([pf_center, pf_center], [0, ymax], 'Color', 'c'); %center of place field
hold off;


% PLOT histogram of #spikes vs. position, regular and normalized (by occupancy time) plots overlaid
axes('position', [0.1 0.635 0.8 0.08]);
if horseshoe_flag == 1
    bincenters = -358:4:358;
else
    bincenters = 1:2:359;
end %if horseshoe
[nn] = hist(V_at_cellts,bincenters);
bar(bincenters, nn, 'k'); hold on;   %spike histo

[occupancy] = hist(Vdata,bincenters);   % occupancy is proportional to number of video frames in bin

if max(nn) ~= 0
    nn_normalized = (nn./occupancy);
    nn_normed_rescaled = nn_normalized * max(nn)/max(nn_normalized);   % rescale so that maxima have same y-value
    plot(bincenters,nn_normed_rescaled,'g'); %firing normalized by occupancy time
    occupancy_rescaled = occupancy*max(nn)/max(occupancy);  %rescale so maxima have same y-value
    ymax = max(nn);
else 
    occupancy_rescaled = occupancy/max(occupancy);  %if no firing data, just scale occupancy relative to itself
    ymax = 1;
end %if max(nn)
plot(bincenters,occupancy_rescaled,'m'); %occupancy time

%line([0,0], [0,ymax], 'Color','r'); %reward site

line([pf_center, pf_center], [0, ymax], 'Color', 'c'); %place field center

if horseshoe_flag == 1
    axis([-360 360 0 ymax]);
else
    axis([0 360 0 ymax]);
end %if horseshoe
%xlabel('Position (degrees)');
ylabel('Spikes');
set(gca, 'XTickLabel', []);
hold off;


% PLOT raster of events, spikes, and occupancy-normalized firing in each lap
height = 1;
bar_fraction = 0.8; 

% find global maximum of Spotsize
% use a 95% percentile as the maximum spotsize, to avoid domination by a few outliers.
ss = [];
for i=1:nLaps
   ix_nzero = find(SpotSize{i} > 0);
   ss = [ss SpotSize{i}(ix_nzero)];     % accumulate nonzero spotsizes
end % for i

switch length(ss)
case 0
    SpotSizeMax = 0;
case 1
    SpotSizeMax = ss(1);
otherwise
    sortss = sort(ss);
    x = [floor(.95*length(ss)) floor(.95*length(ss))+1];
    y = sortss(x);
    SpotSizeMax = interp1(x,y,.95*length(ss));
end % switch

axes('position', [0.1 0.06 0.8 0.568]);
for i = 1:length(SV)
    
  %plot spikes in each lap  
  sp = SV{i};
  sx = [sp sp repmat(NaN, length(sp), 1)];
  sy = repmat([(i*height) (i*height + height *bar_fraction) NaN], length(sp), 1);
  sx = reshape(sx', 1, length(sp)*3);
  sy = reshape(sy', 1, length(sp)*3);
  line(sx, sy, 'Color', [0.7 0.7 0.7]);  hold on;
  
  %plot events in each lap
  stx1 = eventV1{i};
  sty1 = ones(size(stx1))*i +0.5;
  stx2 = eventV2{i};
  sty2 = ones(size(stx2))*i +0.5;
  plot(stx1, sty1,'rx','MarkerSize',6,'LineWidth',1.5 ); 
  plot(stx2, sty2,'gx','MarkerSize',6,'LineWidth',1.5  );
  
  %plot "spots" (occupancy-normalized firing) in each lap
  SpotMaxPlotSize = 7; 
  %%SpSz = SpotMaxPlotSize*SpotSize{i}./(max(SpotSize{i})+1/10000);  % normalize by maximum spotsize of each lap
  SpSz = SpotMaxPlotSize*SpotSize{i}./(SpotSizeMax+1/10000);         % normalize by global maximum 
  SpSz(find(SpSz > SpotMaxPlotSize)) = SpotMaxPlotSize;         % cut off outliers at 15*SpotSizeMax
  spotsY = ones(size(occcenters))*i +0.5;
  noplot = [];
  for j=1:length(occcenters)
      if occcenters(j)<=Vnorun_ends{i}(2) | occcenters(j)>=Vnorun_ends{i}(4)
          noplot = 1;
      end %if 
      if horseshoe_flag==1
        if ~isempty(Vnorun_mid{i}) 
            if occcenters(j)>=Vnorun_mid{i}(1) & occcenters(j)<=Vnorun_mid{i}(2)
                noplot = 1;
            end %if occ...
        end %if ~isempty
      end %if horseshoe
      if isempty(noplot)    
        plot(occcenters(j),spotsY(j),'ok','MarkerSize',max(SpSz(j),0.1))
      end %if isempty
      noplot = [];
  end % for j
  
  %plot magenta horizontal lines where no running occurred
  y_norun_ends = [i+.5 i+.5 NaN i+.5 i+.5];
  line(Vnorun_ends{i},y_norun_ends,'Color','m','LineWidth',1);
  if horseshoe_flag == 1
    if ~isempty(Vnorun_mid{i})
        y_norun_mid = [i+.5 i+.5];
        line(Vnorun_mid{i},y_norun_mid,'Color','m','LineWidth',1);
    end %if ~isempty
  end %if horseshoe

end % for i


line([pf_center, pf_center], [0, nLaps+1], 'Color', 'c'); %place field center

if horseshoe_flag == 1
    axis([-360 360 0 nLaps+1]);
else
    axis([0 360 0 nLaps+1]);
end %if horseshoe
xlabel('Position (degrees)');
ylabel('Lap');
if length(varargin) > 0
    set(gca, 'YTick', 1:nLaps, 'YTickLabel', lapnums);
end %if
hold off;
orient tall;
