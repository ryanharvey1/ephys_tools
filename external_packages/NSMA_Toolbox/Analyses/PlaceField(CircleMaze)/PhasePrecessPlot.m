function fh = PhasePrecessPlot(V,PHtsd,eeg_theta,horseshoe_flag,pf_center,plotrange,eeg_scalefactor,titlestr)

% PhasePrecessPlot  Produces phase precession plot- theta phase of spikes vs. position of animal relative to place field
%
% fh = PhasePrecessPlot(V,PHtsd,eeg_theta,horseshoe_flag,pf_center,plotrange,eeg_scalefactor,titlestr)
%
% INPUTS:
%   V - position tsd; angle on circular track in degrees [0..360) 
%   PHtsd - phase tsd of one cell; data  is a 2 column array with [phase, cyclenum] pairs. Phase is in range 0..1 units 
%   eeg_theta - tsd of eeg filered for theta
%   horseshoe_flag - 1=animal runs in "horseshoe" pattern on circular track, where it runs once
%       around circle in one direction, turns and runs once around in reverse direction, etc.
%       0(or no inputs)=animal runs regular laps around circle
%   pf_center - location (degrees) of center of place field
%   plotrange - x-axis range [-plotrange, plotrange] in theta cycle numbers around placefield center at 0
%   eeg_scalefactor - scaling factor to convert Cheetah A2D values [-2048,2047] to microvolts:
%       Formula: eeg_scalefactor [microvolts/A2Dstep] = 10^6 * 4 [Volts]/(4096*AmplifierGainFactor*A2DGainFactor)
%       AmplifierGainFactor = Vaule set in cheetah control panel or config file
%       A2DgainFactor = A2D programmable gain factor setting (in steps of 1,2,4,8)
%   titlestr - tetrode & cell no. to be used in fig title
% OUTPUTS:
%   fh - figure handle
%
% Theta phase plotted 2x on y-axis in order to see "banana-shaped" pattern of precession better
%
% PL '01, last modified '04 by MN


Vdata = Data(V);
Vts = Range(V,'ts');

if ~isempty(data(PHtsd))
    
    % Reparameterize to degrees in horseshoe pattern (-360->0->360) and count laps.  
    % Check for time gaps in position data - if > 20 sec, make a break there and when data resumes
    % after gap, count as a new lap
    Vts = Range(V,'ts');
    ixbeginepochs = [1; find(diff(Vts) > 20*10000)+1];
    tsbeginepochs = Vts(ixbeginepochs);
    VtsTemp = [];  Vdata = [];  flagixlaps = [];  lastlapinepochs = [];  cumlaps = 0;
    for i = 1:length(ixbeginepochs)
        if i<length(ixbeginepochs)
            Vepoch = Restrict(V,Vts(ixbeginepochs(i)),Vts(ixbeginepochs(i+1)-1));
        else
            Vepoch = Restrict(V,Vts(ixbeginepochs(i)),Vts(end));
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
    
    
    % Limit everything to just running epochs (M1,M2,etc., but not S).  
    % Restrict to full laps
    ixlastlaps = ix_laps(lastlapinepochs);
    eeg_theta = Restrict(eeg_theta,tsbeginepochs,Vts(ixlastlaps));
    PHtsd = Restrict(PHtsd,tsbeginepochs,Vts(ixlastlaps));
    PH = data(PHtsd);
    cell_ts = Range(PHtsd,'ts');
    flagixlaps = data(Restrict(tsd(Vts,flagixlaps),tsbeginepochs,Vts(ixlastlaps)));
    ix_laps = find(flagixlaps == 1);
    V = Restrict(V,tsbeginepochs,Vts(ixlastlaps));
    Vts = Range(V,'ts');  Vdata = data(V);
    V_at_cellts = data(Restrict(V,cell_ts));


    % find timestamp of position closest to the center of place field at each lap
    ixLast = 0;
    ts_center = zeros(nLaps,1);
    for iL = 1:nLaps
        ts1 = Vts(ixLast+1);
        ts2 = Vts(ix_laps(iL));
        [Vmin, ix_center] = min( abs( Vdata(ixLast+1:ix_laps(iL))-pf_center ) );
        ts_center(iL) = Vts(ixLast+ix_center);
        ixLast = ix_laps(iL);    
    end % for iL


    %find indices of first spikes in each cycle
    ix_firstSpikes = 1+find( diff(PH(:,2)) > 0);
    

    %main fig window
    sc = [0 0 1 1];
    fh = figure;
    set(fh,'Units','normalized','Position',[sc(1)+.05,sc(2)+.1,sc(3)*0.9,sc(4)*0.8]);
    orient tall;

    % plot lap-averaged theta-filtered EEG
    timerange = 5;                  % theta time halfwindow in sec
    eegts = Range(eeg_theta,'ts');
    eeg = Data(eeg_theta);
    dT = mean(diff(eegts))/10000;   % delta T between eeg timestamps in sec
    Tbins = -timerange:dT:timerange;                % range of eeg time window
    nTbins = length(Tbins);
    nTbins = floor(nTbins/2)*2;       % make nTbins even 
    Tbins = Tbins(1:nTbins);         % truncate to even if necessary
    avgEEG = zeros(nTbins,1);
    eeg = [avgEEG; eeg; avgEEG];       % shift eeg by nTbins and add zeros at beginning and end
    j = 0;
    for iL=1:nLaps
        % find the index of theta cycle top closest to the place field center crossing
        ixeeg = find(eegts >= ts_center(iL)-dT*5000)+nTbins;        % ixeeg(1) is the closest eeg timestamp to ts_center{iL}
        if  ~isempty(ixeeg)
            ic = ixeeg(1);  
            nxwin = floor(1/(2*dT));                                    
            ixwin = [(ic-nxwin+1):(ic+nxwin)];                % indices of a [-0.5,0.5] sec window around ic
            dfweegwin = diff([0; diff([0; eeg(ixwin)]) ]);    % 2nd forward derivative 
            dbweegwin = [dfweegwin(2:end); 0];                % shift to left by 1 index position
            ix_thetatops = find(dfweegwin >0 & dbweegwin <0); % indices of theta maximums are zero crossings of 2nd derivative
            [thtop, iclosest_theta_top] = min(abs(ix_thetatops - nxwin)); 
            ithtop = ic  + (ix_thetatops(iclosest_theta_top) - nxwin);
            avgEEG = avgEEG + eeg(ithtop-nTbins/2+1:ithtop+nTbins/2);
        else
            j = j+1;
        end %if
    end % for iL
    avgEEG = eeg_scalefactor*avgEEG/(nLaps-j);     % normalize to mean

    axes('position', [0.1 0.885 0.8 0.07]);
    plot(Tbins, avgEEG);
    
    y0 = min(avgEEG);
    y1 = max(avgEEG);
    
    % plot vertical line indicating center of place field
    line([0, 0], [y0, y1], 'Color','c');
    
    axis([-timerange timerange y0 y1]);
    ylabel('avg. \theta-EEG');
    title(titlestr);
    %text(3.5,y0+(y1-y0)/5,'time (sec)');
    xlabel('time (sec)');
    hold off;

    % overall histogram of #spikes vs. position, regular and normalized (by occupancy time) plots overlaid
    axes('position', [0.1 0.735 0.8 0.07]);
    if horseshoe_flag == 1
        bincenters = -358:4:358;
    else
        bincenters = 1:2:359;
    end %if horseshoe
    [nn] = hist(V_at_cellts,bincenters);
    bar(bincenters, nn); hold on;   %spike histo
    
    [occupancy] = hist(Vdata,bincenters);   % occupancy is proportional to number of video frames in bin
    
    nn_normalized = (nn./occupancy);
    nn_normed_rescaled = nn_normalized * max(nn)/max(nn_normalized);   % rescale so that maxima have same y-value
    plot(bincenters,nn_normed_rescaled,'g'); %firing normalized by occupancy time
    
    occupancy_rescaled = occupancy*max(nn)/max(occupancy);  %rescale so maxima have same y-value
    plot(bincenters,occupancy_rescaled,'m'); %occupancy time

    %line([0,0], [0,max(nn)], 'Color','r'); %reward site
    
    % plot vertical line indicating center of place field
    line([pf_center, pf_center], [0, max(nn)], 'Color','c');

    if horseshoe_flag == 1
        axis([-360 360 0 max(nn)]);
    else
        axis([0 360 0 max(nn)]);
    end %if horseshoe
    set(gca, 'XTickLabel', []);
    ylabel('# Spikes');
    hold off;


    %PLOT phase vs postition plots
    %axes('position',[0.1 0.4 0.8 0.32]);
    axes('position',[0.1 0.44 0.8 0.29]);
    plot([V_at_cellts; V_at_cellts], [PH(:,1); 1 + PH(:,1)] ,'.','MarkerSize',4);  %original markersize .2
    hold on;
    plot([V_at_cellts(ix_firstSpikes); V_at_cellts(ix_firstSpikes)], [PH(ix_firstSpikes,1); 1 + PH(ix_firstSpikes,1)] ,'.r','MarkerSize',4);
    % plot vertical line indicating center of place field
    line([pf_center, pf_center], [0, 2], 'Color','c');
    
    if horseshoe_flag == 1
        axis([-360,360,0,2]);
    else
        axis([0,360,0,2]);
    end %if horseshoe
    xlabel('Position (degrees)');
    ylabel('Phase');


    %zoomed-in version centered at place field
    axes('position',[0.1 0.09 0.8 0.29]);
    plot([V_at_cellts; V_at_cellts], [PH(:,1); 1 + PH(:,1)] ,'.','MarkerSize',4);  %original markersize .2
    hold on;
    plot([V_at_cellts(ix_firstSpikes); V_at_cellts(ix_firstSpikes)], [PH(ix_firstSpikes,1); 1 + PH(ix_firstSpikes,1)] ,'.r','MarkerSize',4);
    % plot vertical line indicating center of place field
    line([pf_center, pf_center], [0, 2], 'Color','c');
    
    axis([pf_center-plotrange,pf_center+plotrange,0,2]);
    xlabel('Position (degrees - zoomed in around place field)');
    ylabel('Phase');
    
    
else % if data(PHtsd) is empty (no spikes during maze)
    
    sc = [0 0 1 1];
    fh = figure;
    set(fh,'Units','normalized','Position',[sc(1)+.05,sc(2)+.1,sc(3)*0.9,sc(4)*0.8]);
    orient tall;
    title(titlestr);
    text(.3,.5,'(No Spikes During Maze)','Color','r','FontSize',20);
    set(gca,'XTick',[]); set(gca,'XTickLabel',{});
    set(gca,'YTick',[]); set(gca,'YTickLabel',{});
    

end % if ~isempty(data(PHtsd))




