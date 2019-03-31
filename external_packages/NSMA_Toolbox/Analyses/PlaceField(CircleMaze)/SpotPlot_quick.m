function SpotPlot_quick(V, cell_ts, titlestr, horseshoe_flag, event_ts1, event_ts2)

% SpotPlot_quick  Produces "SpotPlot" for user of 'Run_PlaceFieldGraphs.m' to determine location of place field
% 
% SpotPlot_quick(V, cell_ts, titlestr, horseshoe_flag, event_ts1, event_ts2)
% 
% INPUTS:
%   V - position tsd; angle on circular track in degrees [0..360) 
%   cell_ts - array of timestamps of spikes from one cell
%   titlestr - tetrode & cell no. to be used in fig title
%   horseshoe_flag - 1=animal runs in "horseshoe" pattern on circular track, where it runs once
%       around circle in one direction, turns and runs once around in reverse direction, etc.
%       0(or no inputs)=animal runs regular laps around circle
%   event_ts1 & 2 - array of timestamps of events for up to 2 different types of events.  These
%       inputs optional - no events need be entered, or inputting just one type of event (_ts1)
%       but not _ts2 is also acceptible.
% OUTPUTS:
%   (none)
%
% PL '01, last modified '04 by MN


% Count Laps
% Check for time gaps in position data - if > 20 sec, make a break there and when data resumes
% after gap, count as a new lap
% Vts = Range(V,'ts');
Vts = range(V,'includenan'); %Ryan Harvey 12/9/16
ixbeginepochs = [1; find(diff(Vts) > 20*10000)+1];
tsbeginepochs = Vts(ixbeginepochs);
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


% Limit cell_ts and event_ts to just running epochs (M1,M2,etc., but not S1,S2).  
% Restrict to full laps - don't use the last partial lap
ixlastlaps = ix_laps(lastlapinepochs);
cell_ts = data(restrict(ts(cell_ts),tsbeginepochs,Vts(ixlastlaps)));
event_ts1 = data(restrict(ts(event_ts1),tsbeginepochs,Vts(ixlastlaps)));
event_ts2 = data(restrict(ts(event_ts2),tsbeginepochs,Vts(ixlastlaps)));
flagixlaps = data(restrict(tsd(Vts,flagixlaps),tsbeginepochs,Vts(ixlastlaps)));
ix_laps = find(flagixlaps == 1);
V = restrict(V,tsbeginepochs,Vts(ixlastlaps));
Vts = Range(V,'ts');  Vdata = data(V);


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
    ixLast = ix_laps(iL);
end % for iL


sc = [0 0 1 1];
fh = figure;
set(fh,'Units','normalized','Position',[sc(1)+.05,sc(2)+.4,sc(3)*0.9,sc(4)*0.5]);
orient tall;


% PLOT histogram of #spikes vs. position, regular and normalized (by occupancy time) plots overlaid
axes('position', [0.1 0.81 0.8 0.1]);
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

if horseshoe_flag == 1
    axis([-360 360 0 ymax]);
else
    axis([0 360 0 ymax]);
end % if horseshoe

%xlabel('Position (degrees)');
ylabel('Spikes');
title(titlestr);
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

axes('position', [0.1 0.1 0.8 0.7]);
for i = 1:length(SV)
  sp = SV{i};
  sx = [sp sp repmat(NaN, length(sp), 1)];
  sy = repmat([(i*height) (i*height + height *bar_fraction) NaN], length(sp), 1);
  sx = reshape(sx', 1, length(sp)*3);
  sy = reshape(sy', 1, length(sp)*3);
  
  stx1 = eventV1{i}; 
  sty1 = ones(size(stx1))*i +0.5;
  stx2 = eventV2{i};
  sty2 = ones(size(stx2))*i +0.5;
  line(sx, sy, 'Color', [0.7 0.7 0.7]);
  hold on;
  plot(stx1, sty1,'rx','MarkerSize',6,'LineWidth',1.5 ); 
  plot(stx2, sty2,'gx','MarkerSize',6,'LineWidth',1.5  );
  SpotMaxPlotSize = 7; 
  %%SpSz = SpotMaxPlotSize*SpotSize{i}./(max(SpotSize{i})+1/10000);  % normalize by maximum spotsize of each lap
  SpSz = SpotMaxPlotSize*SpotSize{i}./(SpotSizeMax+1/10000);         % normalize by global maximum 
  SpSz(find(SpSz > SpotMaxPlotSize)) = SpotMaxPlotSize;         % cut off outliers at 15*SpotSizeMax
  spotsY = ones(size(occcenters))*i +0.5;
  for j=1:length(occcenters)
      plot(occcenters(j),spotsY(j),'ok','MarkerSize',max(SpSz(j),0.1))
  end % for j
end % for i

if horseshoe_flag == 1
    axis([-360 360 0 nLaps+1]);
else
    axis([0 360 0 nLaps+1]);
end %if horseshoe

xlabel('Position (degrees)');
ylabel('Lap');
hold off;
orient tall;
