function fh = PhaseTimePlot(V,PHtsd,eeg_theta,horseshoe_flag,event_ts1,event_ts2,pf_center,plotrange,titlestr)

% PhaseTimePlot  Produces lap-by-lap raster of "theta cycle time" of spikes occuring around place field
%
% fh = PhaseTimePlot(V,PHtsd,eeg_theta,horseshoe_flag,event_ts1,event_ts2,pf_center,plotrange,titlestr)
%
% INPUTS:
%   V - position tsd; angle on circular track in degrees [0..360)
%   PHtsd - phase tsd of one cell; data  is a 2 column array with [phase, cyclenum] pairs. Phase is in range 0..1 units
%   eeg_theta - tsd of eeg filered for theta
%   horseshoe_flag - 1=animal runs in "horseshoe" pattern on circular track, where it runs once
%       around circle in one direction, turns and runs once around in reverse direction, etc.
%       0(or no inputs)=animal runs regular laps around circle
%   event_ts1 & 2 - array of timestamps of events for up to 2 different types of events.  These
%       inputs optional (may be empty [])
%   pf_center - location (cm) of center of place field
%   plotrange - x-axis range [-plotrange, plotrange] in theta cycle numbers around placefield center at 0
%   titlestr - tetrode & cell no. to be used in fig title
% OUTPUTS:
%   fh - figure handle
%
% Lap numbers on y-axis;  On x-axis: theta cycles before and after animal passes through center of place field
%   (center of place field: theta cycle = 0)
%
% PL '01, last modified '04 by MN
%
% New modifications to handle multiple fields by Drew Maurer
% I've rigged this code for the Precession\Assemblies 2006 paper so that
% the variables expected are:
% INPUTS:
% V- position tsd; angle (IN DEGREES) on circular track [0..360]
% PHtsd - phase tsd of one cell; data  is a 2 column array with [phase, cyclenum] pairs. Phase is in range 0..1 units
%             *************This input IS unecessary! Instead one can use a
%             list of timestamps. Comment out line 56 and uncomment line 60
% eeg_theta - tsd of eeg filered for theta. This variable can be generated
%             by Filter4Theta.
% horseshoe_flag - 1=animal runs in "horseshoe" pattern on circular track, where it runs once
%       around circle in one direction, turns and runs once around in reverse direction, etc.
%       0(or no inputs)=animal runs regular laps around circle.
%       NOW IT WORKS. V_epoch was created by the LapCount code, but wasn't
%       being introduced into the code.
% event_ts1 & 2 - array of timestamps of events for up to 2 different 
%       types of events. I've been using the spike timestamps from the other two place fields. Now they are plotted as lines instead of circles or * 
%       inputs are optional (may be empty []).
% pf_center - THIS CODE DOES NOT ACCEPT PF CENTER IN CM!!!!!! Instead
%       Miriam wanted the pf_center in degrees! Convert you pf_center(cm) to
%       degrees by: pf_center = (max(Data(V)).*pf_center)./(INSERT THE CIRCUMFERENCE OF YOUR TRACK HERE);
% plotrange - x-axis range [-plotrange, plotrange] in theta cycle
%       numbers; this code expects two numbers [Numbers of cycles prior to
%       center (negative), Number of cycles following center]. I've been using
%       [-15 15]


my_colors = {'r' 'b' 'g' 'm' 'k' 'y' 'c' 'r' 'b' 'g' 'm' 'k' 'y' 'c' 'r' 'b' 'g' 'm' 'k' 'y' 'c'  };

PH= data(PHtsd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is how one could put in a list of timestamps instead of the code
% requiring you provide the spike phase:
%    PH = ThetaPhase({ts(PHtsd)},eeg_theta,PHtsd(1,1) -20000 ,PHtsd(end)+20000); PH= data(PH);
%    Where PHtsd = the first required input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(PH)

    % Reparameterize to degrees in horseshoe pattern (-360->0->360) and count laps.
    % Check for time gaps in position data - if > 20 sec, make a break there and when data resumes
    % after gap, count as a new lap
    Vts = Range(V,'ts');
    % This code is asking for epoch indices of laps
    %     ixbeginepochs = [1; find(diff(Vts) > 20*10000)+1];
    %     tsbeginepochs = Vts(ixbeginepochs);
    %     ixendepochs = [find(diff(Vts) > 20*10000); length(Vts)];
    %     tsendepochs = Vts(ixendepochs);
    if horseshoe_flag == 1;
        %load([pwd '\spike_rast_vars_big.mat']);
        ixbeginepochs = find(Vts(1) == Vts);
        tsbeginepochs = Vts(ixbeginepochs);
        ixendepochs = find(Vts(end) == Vts);
        tsendepochs = Vts(ixendepochs);
    end
    if horseshoe_flag == 0;
        % load([pwd '\spike_rast_vars_small.mat']);
        ixbeginepochs = find(Vts(1) == Vts);
        tsbeginepochs = Vts(ixbeginepochs);
        ixendepochs = find(Vts(end) == Vts);
        tsendepochs = Vts(ixendepochs);
    end
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
        VtsTemp = [VtsTemp; Range(Vepoch,'ts')];  Vdata = [Vdata; Data(V_epoch)];
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
    PHtsd = restrict(PHtsd,tsbeginepochs,ts_end);
    PH = data(PHtsd);
    cell_ts = Range(PHtsd,'ts');
    flagixlaps = data(restrict(tsd(Vts,flagixlaps),tsbeginepochs,ts_end));
    ix_laps = find(flagixlaps == 1); nLaps = length(ix_laps);
    V = restrict(V,tsbeginepochs,ts_end);
    Vts = Range(V,'ts');  Vdata = data(V);
    PH_at_cellts = PH(:,2)+PH(:,1);     % theta time of spikes
    PH_events1 = ThetaPhase({ts(event_ts1)},eeg_theta,tsbeginepochs,ts_end);
    event_ts1 = Range(PH_events1{1},'ts');
    PH_events1_data = data(PH_events1{1});
    if ~isempty(PH_events1_data)
        PH_at_eventts1 = PH_events1_data(:,2)+PH_events1_data(:,1);    % theta time of events
    else
        PH_at_eventts1 = [];
    end %if
    PH_events2 = ThetaPhase({ts(event_ts2)},eeg_theta,tsbeginepochs,ts_end);
    event_ts2 = Range(PH_events2{1},'ts');
    PH_events2_data = data(PH_events2{1});
    if ~isempty(PH_events2_data)
        PH_at_eventts2 = PH_events2_data(:,2)+PH_events2_data(:,1);    % theta time of events
    else
        PH_at_eventts2 = [];
    end %if
    PH_endepochs = ThetaPhase({ts(ts_end)},eeg_theta,tsbeginepochs,ts_end);
    PH_endepochs_data = data(PH_endepochs{1});
    PH_at_endepochs = PH_endepochs_data(:,2)+PH_endepochs_data(:,1);


    % find theta cycle number closest to the center of place field at each lap
    ixLast = 0;
    ts_center = zeros(nLaps,1);
    norun_laps = [];
    for iL = 1:nLaps
        ts1 = Vts(ixLast+1);
        ts2 = Vts(ix_laps(iL));
        [Vmin, ix_center] = min( abs( Vdata(ixLast+1:ix_laps(iL))-pf_center ) );
        ts_center(iL) = Vts(ixLast+ix_center);
        %record laps where rat ends running epoch before reaching pf
        if pf_center > max(Vdata(ixLast+1:ix_laps(iL)))
            norun_laps = [norun_laps; iL];
        end %if
        ixLast = ix_laps(iL);
    end % for iL


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % I put the code in below to make sure the spikes were being plotted
    % onto position correctly- It looks like they are
    % The ts_center value is not, however!
    %
    % Epilouge: The problem was that the reparameterzed data from the
    % Horseshoe track, "V_epoch" (line 54) wasn't making it into the code.
    % This was fixed by switching Vepoch with V_epoch at line 59.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     cell_pos = interp1(Vts, Vdata,cell_ts);
    %     figure; plot(Vts, Vdata)
    %     hold on; plot(cell_ts,cell_pos, '.')
    %     hold on; plot(Vts(ix_laps), Vdata(ix_laps), 'r.')
    %     hold on; plot(ts_center, 225, 'm.')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PH_center = ThetaPhase({ts(ts_center)},eeg_theta,tsbeginepochs,ts_end);
    PHdata_center = data(PH_center{1});
    cycl_center = round(PHdata_center(:,2) + PHdata_center(:,1));


    % caluclate Phase times per Lap
    ixLast = 0;
    all_cell_phase_times = [];
    all_event1_phase_times = [];
    all_event2_phase_times = [];
    for iL = 1:nLaps
        reach_pf = 1;
%         if ~isempty(norun_laps)
%             if ~isempty(find(norun_laps == i))
%                 reach_pf = 0;
%             end %if
%         end %if
        if reach_pf %if rat reaches pf in that lap
            ts1 = Vts(ixLast+1)-60*10000;
            ts2 = Vts(ix_laps(iL))+60*10000;
            cell_phase_times{iL} = PH_at_cellts(find(cell_ts >= ts1 & cell_ts <= ts2) ) - cycl_center(iL);
            all_cell_phase_times = [all_cell_phase_times; cell_phase_times{iL}];   % append cell phases for current lap
            event1_phase_times{iL} = PH_at_eventts1(find(event_ts1 >= ts1 & event_ts1 <= ts2) ) - cycl_center(iL);
            all_event1_phase_times = [all_event1_phase_times; event1_phase_times{iL}]; % append Event1 phases for current lap
            event2_phase_times{iL} = PH_at_eventts2(find(event_ts2 >= ts1 & event_ts2 <= ts2) ) - cycl_center(iL);
            all_event2_phase_times = [all_event2_phase_times; event2_phase_times{iL}]; % append Event2 phases for current lap
            endepoch_phase_times{iL} = PH_at_endepochs(find(ts_end >= ts1 & ts_end <= ts2) ) - cycl_center(iL);
        end %if
        ixLast = ix_laps(iL);
    end % for iL


    % overall histogram of #spikes vs. theta time
    bincenters = plotrange(1):0.1:plotrange(2);
    [nn] = hist(all_cell_phase_times,bincenters);
    [xx] = hist(all_event1_phase_times,bincenters);
    [yy] = hist(all_event2_phase_times,bincenters);

    nn([1 end]) = [];                 % remove overflow bins
    xx([1 end]) = [];                 % remove overflow bins
    yy([1 end]) = [];                 % remove overflow bins
    bincenters([1 end]) = [];

    sc = [0 0 1 1];
    fh = figure;
    set(fh,'Units','normalized','Position',[sc(1)+.05,sc(2)+.1,sc(3)*0.9,sc(4)*0.8]);
    orient tall;

    if length(event_ts2) > 1;
        subplotboxes = 12;
    else
        subplotboxes = 10;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %axes('position', [0.1 0.8 0.8 0.1]);
    subplot(subplotboxes, 1, 1)
    bar(bincenters, nn, my_colors{1});

    if max(nn) > 0
        ymin = -max(nn)/4;
        ymax = 1.25*max(nn);
    else
        ymin = -.25;
        ymax = 1.25;
    end % if max(nn)

    hold on;
    % plot vertical line indicating center of place field
    xline = plotrange(1):1:plotrange(2);
    xline = [xline;xline];
    yline = [ymin, ymax]
    Yline = repmat(yline', 1, length(xline))
    plot(xline,Yline, 'k')
    line([0, 0], [ymin, ymax], 'Color','c');

    xrange = plotrange(1):plotrange(2);
    set(gca,'xtick',xrange);
    %for i=1:5:length(xrange)
    %    cellarray{i} = xrange(i);
    %end
    %set(gca,'XTickLabel',cellarray);
    set(gca, 'XTickLabel', []);

    ylabel('# Spikes');
    %xlabel('Phase-Time (CycleNum.Phase)');
    axis([plotrange(1) plotrange(2) ymin ymax]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold on;
    subplot(subplotboxes, 1, 2)
    bar(bincenters, xx, my_colors{2});

    if max(nn) > 0
        ymin = -max(xx)/4;
        ymax = 1.25*max(xx);
    else
        ymin = -.25;
        ymax = 1.25;
    end % if max(nn)
    hold on;
    % plot vertical line indicating center of place field
    xline = plotrange(1):1:plotrange(2);
    xline = [xline;xline];
    yline = [ymin, ymax]
    Yline = repmat(yline', 1, length(xline))
    plot(xline,Yline, 'k')
    line([0, 0], [ymin, ymax], 'Color','c');
    axis([plotrange(1) plotrange(2) ymin ymax]);

    xrange = plotrange(1):plotrange(2);
    set(gca,'xtick',xrange);
    %for i=1:5:length(xrange)
    %    cellarray{i} = xrange(i);
    %end
    %set(gca,'XTickLabel',cellarray);
    set(gca, 'XTickLabel', []);

    ylabel('# Spikes');
    %xlabel('Phase-Time (CycleNum.Phase)');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if subplotboxes == 12;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %for the third field
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hold on;
        subplot(subplotboxes, 1, 3)
        bar(bincenters, yy, my_colors{3});

        if max(nn) > 0
            ymin = -max(yy)/4;
            ymax = 1.25*max(yy);
        else
            ymin = -.25;
            ymax = 1.25;
        end % if max(nn)
        hold on;
        % plot vertical line indicating center of place field
        xline = plotrange(1):1:plotrange(2);
        xline = [xline;xline];
        yline = [ymin, ymax]
        Yline = repmat(yline', 1, length(xline))
        plot(xline,Yline, 'k')
        line([0, 0], [ymin, ymax], 'Color','c');
        axis([plotrange(1) plotrange(2) ymin ymax]);

        xrange = plotrange(1):plotrange(2);
        set(gca,'xtick',xrange);
        %for i=1:5:length(xrange)
        %    cellarray{i} = xrange(i);
        %end
        %set(gca,'XTickLabel',cellarray);
        set(gca, 'XTickLabel', []);

        ylabel('# Spikes');
        %xlabel('Phase-Time (CycleNum.Phase)');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    title(titlestr);
    hold off;


    % raster plot of events & spikes in each lap (over theta time)
    height = 1;
    bar_fraction = 0.8;

    %axes('position', [0.1 0.09 0.8 0.7]);
    if subplotboxes == 12;
        subplot(subplotboxes, 1, 4:subplotboxes)
    else
        subplot(subplotboxes, 1, 3:subplotboxes)
    end
    for i = 1:nLaps
        reach_pf = 1;
        if ~isempty(norun_laps)
            if ~isempty(find(norun_laps == i))
                reach_pf = 0;
            end %if
        end %if
        if reach_pf %if rat reaches pf in that lap
            sp = cell_phase_times{i};
            sx = [sp sp repmat(NaN, length(sp), 1)];
            sy = repmat([(i*height) (i*height + height *bar_fraction) NaN], length(sp), 1);
            sx = reshape(sx', 1, length(sp)*3);
            sy = reshape(sy', 1, length(sp)*3);

            stx = event1_phase_times{i};
            stx1 = [stx stx repmat(NaN, length(stx), 1)];
            sty1 = repmat([(i*height) (i*height + height *bar_fraction) NaN], length(stx), 1);
            stx1 = reshape(stx1', 1, length(stx)*3);
            sty1 = reshape(sty1', 1, length(stx)*3);

            stxtwo = event2_phase_times{i};
            %sty2 = ones(size(stx2))*i +0.5;
            stx2 = [stxtwo stxtwo repmat(NaN, length(stxtwo), 1)];
            sty2 = repmat([(i*height) (i*height + height *bar_fraction) NaN], length(stxtwo), 1);
            stx2 = reshape(stx2', 1, length(stxtwo)*3);
            sty2 = reshape(sty2', 1, length(stxtwo)*3);

            line(sx, sy, 'Color', my_colors{1});
            hold on;
            plot(stx1, sty1,my_colors{2});
            plot(stx2, sty2,my_colors{3});
            if ~isempty(endepoch_phase_times{i})
                line([endepoch_phase_times{i} plotrange(2)],[i+.5 i+.5],'Color','m','LineWidth',1);
            end %if ~isempty
        else %if rat ends running epoch before reaching pf in that lap
            line([plotrange(1) plotrange(2)],[i+.5 i+.5],'Color','m','LineWidth',1);
        end %if isempty
    end % for i

    % plot vertical line indicating center of place field
    xline = plotrange(1):1:plotrange(2);
    xline = [xline;xline];
    yline = [0, nLaps+1];
    Yline = repmat(yline', 1, length(xline));
    plot(xline,Yline, 'k');
    line([0, 0], [0, nLaps+1], 'Color','c');

    axis([plotrange(1) plotrange(2) 0 nLaps+1]);

    set(gca,'xtick',xrange);
    for i=1:5:length(xrange)
        cellarray{i} = xrange(i);
    end
    set(gca,'XTickLabel',cellarray);

    ylabel('Lap');
    xlabel('Phase-Time (CycleNum.Phase)');
    hold off;
    orient tall;


else  % if PH is empty (no spikes during maze)

    sc = [0 0 1 1];
    fh = figure;
    set(fh,'Units','normalized','Position',[sc(1)+.05,sc(2)+.1,sc(3)*0.9,sc(4)*0.8]);
    orient tall;
    title(titlestr);
    text(.3,.5,'(No Spikes During Maze)','Color','r','FontSize',20);
    set(gca,'XTick',[]); set(gca,'XTickLabel',{});
    set(gca,'YTick',[]); set(gca,'YTickLabel',{})


end % if ~isempty(PH)
