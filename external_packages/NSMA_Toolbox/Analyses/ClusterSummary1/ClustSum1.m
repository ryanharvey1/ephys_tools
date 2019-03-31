function [f] = ClustSum1(epochs, t_file, tt_file, varargin)

% ClustSum1  Creates cluster summary ("BPALL") sheet: avg. waveform, ISI hist, autocorr, etc.
%
% [f] = ClustSum1(epochs, t_file, tt_file, varargin)
%
% INPUTS: 
%   epochs - a struct object with fields:
%       epochs.names = cell array of strings of epoch names 
%       epochs.intervals = cell array of 1x2 arrays with [start_ts  end_ts] (start
%           and end timestamps of each epoch) -- elements in cell array correspond
%           to epochs, listed in same sequence as in epochs.names
%   t_file - the name of the tfile. If the name is partial(i.e. TT4) then
%                 check cluster from file will do a search for all files with that
%                 prefix.
%   tt_file - the name of the tt (or lacking that, the dat) file
%   varargin PARAMETERS: 
%       These optional additional parameters allow additional plots to be 
%       added to the check cluster plot. The format for these plots
%       is 'parameter name' (in single quote), parameter values.
%       
%       To plot place field and scatter field plots:
%
%       f = Check_cluster_from_file2('TT10','C:\Data\Stress\6989_05\TT10.dat',...
%           'position','C:\Data\Stress\6989_05\VT1.ascii',...
%           [28261548, 41298161;63768362, 78104458],'save_figures','eps')
%         
%       the first arguement after 'position' is the name of the ascii position file
%       the second is a matrix of start and end times for the epoch.
% 
%       to save files instead of keeping up a figure window, pass in 
%       'save_figures' as parameter 1 and 'eps' or 'fig' or 'bmp' as parameter 2.
%       parameter 2 specifies the file format (see saveas for other options).
%
% OUTPUTS: 
%   f - file handle of figure
%
% Cowen, from ADR CheckCluster, last modified '03 by MN


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE DECLARATIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
draw_position = 0;
draw_PETH = 0;
save_figures = 0;
text_font_size = 7;
figure_count = 1;
n_place_bins = 14; % Default for place fields.
f = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find all files that correspond the the passed in TT file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get rid of the extension if there is any. We'll put it back later.
[tfile_path, name, ext] = fileparts(t_file);
dirinfo = dir([name '*.t']);
if isempty(dirinfo)
    disp('NO .t FILES BY THAT NAME!')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the varargins if there are any-- and draw stuff.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iV = 1;
while iV <= length(varargin)
    
    if size(varargin{iV},2) == 1
        varargin{iV} = varargin{iV}';
    end
    switch varargin{iV}
    case 'position'
        iV = iV + 1;
        draw_position = 1;
        disp('Loading Position File.')     
        n_place_bins = 14;
        position = load(varargin{iV}); 
        iV = iV + 1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Try a 'smart' way of determining if this is 
        %     Cheetah NT or sun position info:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if position(1,1) > 200000
            time_divisor = 1;
        else
            time_divisor = 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Having a position of 0,0 is typically an error in tracking so get rid of it.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        position = position(find(position(:,2)~=0),:);
        position = position(find(position(:,3)~=0),:);
        epoch_times = varargin{iV};
        iV = iV + 1;
        n_epochs = size(epoch_times,1);
        for epoch = 1:n_epochs
            x_pos{epoch} = Restrict(tsd(floor(position(:,1)/time_divisor),position(:,2)),epoch_times(epoch,1),epoch_times(epoch,2));
            y_pos{epoch} = Restrict(tsd(floor(position(:,1)/time_divisor),position(:,3)),epoch_times(epoch,1),epoch_times(epoch,2));
            x_pos{epoch} = Smooth_tsd(x_pos{epoch},30);
            y_pos{epoch} = Smooth_tsd(y_pos{epoch},30);
        end
    case 'n_place_bins'
        iV = iV + 1;
        n_place_bins = varagin{iV};
        iV = iV + 1;
        
    case 'PETH'
        iV = iV + 1;
        draw_peth = 1;
    case 'save_figures'
        iV = iV + 1;
        save_figures = 1;
        fig_string = varargin{iV};
        iV = iV + 1;
    otherwise
        error('Unknown Parameter')
    end
end



for current_t_file = 1:length(dirinfo)
    current_t_file_name = [tfile_path filesep dirinfo(current_t_file).name];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    spike_times = LoadSpikes({current_t_file_name});
    spike_times = spike_times{1};
    [t, wv] = ReadTT(tt_file,Data(spike_times),1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If no spikes, go to the next file.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(t)
        break
    end
    t = floor(t);
    ix = find( (floor(Data(spike_times))-t) ~= 0);
    if ~isempty(ix)
        warning(['WARNING: ' num2str(length(ix)) ' spikes in the t file do not correspond to spikes in the tt file. Spikes are removed!' ])
        t(ix) = [];
        wv(ix,:,:) = [];
    end
    
    [peak, ipeak] = max(wv, [], 3);
    [vlly, ivlly] = min(wv, [], 3);
    
    
    [junk, tfile_name] = fileparts(current_t_file_name);
    [junk, ttfile_name] = fileparts(tt_file);

    f(figure_count) = figure('Name', tfile_name, 'NumberTitle', 'Off');
    sc = [0 0 1 1];
    set(f(figure_count),'Units','normalized','Position',[sc(1)+.05,sc(2)+.1,sc(3)*0.9,sc(4)*0.8]);
    orient tall;
    
    figure_count = figure_count + 1;
    
    clf
    colormap(1-gray);
    %nPlot = 4;
    nPlot = 5;
    mPlot = 2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT: AvgWaveform
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(nPlot,mPlot,1);
    CO = get(gca,'ColorOrder');
    for it = 1:4
        mWV = squeeze(mean(wv(:,it,:),1));
        sWV = squeeze(std(wv(:,it,:),1));
        xrange = (34 * (it-1)) + (1:32); 
        hold on;
        h = plot(xrange, mWV);
        set(h, 'Color',CO(it,:))
        h=errorbar(xrange,mWV,sWV); 
        set(h, 'Color',CO(it,:))
    end
    axis off
    axis([0 140 -2100 2100])
    set(gca,'FontSize',text_font_size)
    title('Average Waveform');
    
    hold off
    drawnow

    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT: ISIStats (text)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    subplot(nPlot, mPlot, 2);
    
    nmsgs = 1;
    
    msgstr = {tfile_name};
    
    
    nmsgs = nmsgs+1;
    if isempty(getenv('USER'))
        msgstr{nmsgs} = sprintf('Cut on %s', date);
    else   
        msgstr{nmsgs} = sprintf('Cut by %s on %s', getenv('USER'), date);
    end
    nmsgs = nmsgs+1;
    nSpikes = length(t);
    msgstr{nmsgs} = sprintf('%d spikes ', nSpikes);
    nmsgs = nmsgs+1;
    mISI = mean(diff(t));
    fr = 10000 * nSpikes/(t(end) - t(1));
    
    
    msgstr{nmsgs} = sprintf('firing rate = %.4f spikes/sec ', fr);
    
    
    nmsgs = nmsgs+1;
    
    sw = ivlly - ipeak;
    mSW = mean(sw,1);
    vSW = std(sw, 1);
    msgstr{nmsgs}   = sprintf('spikewidth (Ch1,2) = %4.2f +/- %4.2f   %4.2f +/- %4.2f', ...
        mSW(1), vSW(1), mSW(2), vSW(2));
    msgstr{nmsgs+1} = sprintf('spikewidth (Ch3,4) = %4.2f +/- %4.2f   %4.2f +/- %4.2f', ...
        mSW(3), vSW(3), mSW(4), vSW(4));
    nmsgs = nmsgs+2;
    
    p2v = abs(peak)./abs(vlly);
    mp2v = mean(p2v,1);
    vp2v = std(p2v,1);
    
    
    msgstr{nmsgs}   = sprintf('peak/valley (Ch1,2) = %4.2f +/- %4.2f    %4.2f +/- %4.2f', ...
        mp2v(1), vp2v(1), mp2v(2), vp2v(2));
    msgstr{nmsgs+1} = sprintf('peak/valley (Ch3,4) = %4.2f +/- %4.2f    %4.2f +/- %4.2f', ...
        mp2v(3), vp2v(3), mp2v(4), vp2v(4));
    
    
    nmsgs = nmsgs+2;
    h=text(0,0.5, msgstr);
    set(h,'FontSize',text_font_size)
    axis off
    drawnow
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT: HistISI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(nPlot, mPlot, 3);
    HistISI(ts(t));
    set(gca,'FontSize',text_font_size)
    title('histogram of log(ISI)','FontSize',text_font_size)
    ylabel('nSpikes','FontSize',text_font_size);
    xlabel('msec','FontSize',text_font_size)
    axis tight
    drawnow
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Peak plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(nPlot, mPlot, 4);
    plot(peak(:,2),peak(:,1),'.')
    hold on
    plot(peak(:,2),-peak(:,3),'.')
    plot(-peak(:,4),-peak(:,3),'.')
    plot(-peak(:,4),peak(:,1),'.')
    plot(0,0,'+r')
    H = findobj(gca, 'Type', 'Line');
    set(H, 'MarkerSize', 2)
    axis tight
    axis off
    title(['Peak Plot'],'FontSize',text_font_size)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOT: AutoCorr (AutoCorr.m is VERY slow and peter's autocorr.dll C program does not work yet)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(nPlot, mPlot, 5);
    acorr_bin_msec = 4;
    [histvals,x] = AutoCorr(t, acorr_bin_msec, 250); % This is peter's C autocorr
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %acorr(floor(length(acorr)/2+1)) = 0;    % set 0 lag to 0 for scaling
    % Cannot use bar-- there is some matlab bug that prevents it's use in this subplot.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    plot(x,histvals);			% show acorr
    
    xlabel(['msec (' num2str(acorr_bin_msec) 'msec binsize)'],'FontSize',text_font_size)
    ylabel('rate','FontSize',text_font_size)
    h = title('Autocorrelation','FontSize',text_font_size);
    drawnow
    set(gca,'FontSize',text_font_size)
    % PLOT: AutoCorr (AutoCorr.m is VERY slow and peter's autocorr.dll C program does not work yet)
    subplot(nPlot, mPlot, 6);
    
    
    if 1    % draws Stephen's ISI at t vs. ISI at t+1 scatterplot to detect triplets - 
            % if want to use this plot, change to "if 1" - if not, use "if 0"
            
        d = diff(t)/10; % Convert isi to msec.
        H = plot(d(1:end-1),d(2:end),'.');
        set(H, 'MarkerSize', 2)
        set(gca, 'XScale', 'log');
        set(gca, 'YScale', 'log');
        axis([0 1500 0 1500])
        xlabel('msec (t)','FontSize',text_font_size)
        ylabel('msec (t+1)','FontSize',text_font_size)
        set(gca,'FontSize',text_font_size)
    end
    
    
    if 0   % draws original shortened autocorr bar graph -
           % if want to use this plot, change to "if 1" - if not, use "if 0"
           
        acorr_bin_msec = 1;
        [histvals,x] = autocorr(t, acorr_bin_msec, 50); % This is peter's C autocorr
        
        %acorr(floor(length(acorr)/2+1)) = 0;    % set 0 lag to 0 for scaling
        % Cannot use bar-- there is some matlab bug that prevents it's use in this subplot.
        
        bar(x,histvals);			% show acorr
        
        xlabel(['msec (' num2str(acorr_bin_msec) 'msec binsize)'],'FontSize',text_font_size)
        ylabel('rate','FontSize',text_font_size)
        h = title('Autocorrelation','FontSize',text_font_size);
        drawnow
        set(gca,'FontSize',text_font_size)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PLOTS: Energy Variability vs. time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(nPlot,mPlot,7:8);
    CO = get(gca,'ColorOrder');
    nSpikes = length(t);
    energy = zeros(nSpikes,4);
    for it = 1:4
        tmp = squeeze(wv(:,it,:));
        energy(:,it) = sqrt(sum(tmp.^2,2));
    end
    av_energy = mean(energy);
    delta_energy = sum(energy - repmat(av_energy,nSpikes,1),2);
    
    lb = min(delta_energy);
    ub = max(delta_energy);
    delta_scalefactor = max(abs(lb),abs(ub));
    plot(t/10000/60,delta_energy/delta_scalefactor,'.','MarkerSize', 1);   % scale delta_energy into [-1,1] range
    
    y0 = max(delta_energy/delta_scalefactor);
    y1 = min(delta_energy/delta_scalefactor);
    DrawEpochLines(y0,y1,epochs);

    set(gca,'FontSize',text_font_size)
    xlabel('Time (minutes)','FontSize',text_font_size)
    ylabel('\Delta E','FontSize',text_font_size)

    if (epochs.intervals{end}(2)==0)
        disp('warning; no sleep 2');
    else
        axis([epochs.intervals{1}(1)/10000/60 epochs.intervals{end}(2)/10000/60 y1 y0]);
    end
    
    
    subplot(nPlot,mPlot,9:10);
    
    cosphi_energy = (energy * av_energy')./(sqrt(sum(av_energy.^2,2))*norm(energy));
    %cosphi_energy = (energy * av_energy')./(sqrt(sum(av_energy.^2,2))*sqrt(sum(energy.^2,2)));

    lb = min(cosphi_energy);
    ub = max(cosphi_energy);
    cosphi_scalefactor = max(abs(lb),abs(ub));
    plot(t/10000/60+0.02,cosphi_energy/cosphi_scalefactor,'.g','MarkerSize', 1)    %% scale cosphi_energy into [-1,1] range
    
    y0 = max(cosphi_energy/cosphi_scalefactor);
    y1 = min(cosphi_energy/cosphi_scalefactor);
    DrawEpochLines(y0,y1,epochs);
    
    set(gca,'FontSize',text_font_size)
    xlabel('Time (minutes)','FontSize',text_font_size)
    ylabel('\Delta\Phi','FontSize',text_font_size)

    if (epochs.intervals{end}(2)==0)
        disp('warning; no sleep 2');
    else
        axis([epochs.intervals{1}(1)/10000/60 epochs.intervals{end}(2)/10000/60 y1 y0]);
    end
    
    hold off
    drawnow

    
    if 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %subplot(nPlot, mPlot, 7:8);
    
    
    plot(t/10000/60,peak,'.')
   
    set(gca,'FontSize',text_font_size)
    xlabel('Time (minutes)','FontSize',text_font_size)
    ylabel('Peak','FontSize',text_font_size)
    
    H = findobj(gca, 'Type', 'Line');
    
    set(H, 'MarkerSize', 2)
    axis tight
    
    drawnow
    end %if 0 
    
    
    
    if save_figures
        saveas(f(figure_count-1),[tfile_name '_ClustView'],fig_string);
        close(f(figure_count-1))
    end
   
    if draw_position
        f(figure_count) = figure('Name', ['Place Fields ' tfile_name], 'NumberTitle', 'Off');
        figure_count = figure_count + 1;
        subplot_count = 1;
        for epoch = 1:n_epochs
            subplot(n_epochs*2,2,subplot_count)
            subplot_count = subplot_count + 1;
            spikes = Restrict(spike_times,epoch_times(epoch,1), epoch_times(epoch,2));
            [TC,Occ] = TuningCurves({spikes}, x_pos{epoch}, n_place_bins, y_pos{epoch}, n_place_bins);
            NTC = TC./Occ;
            fq = length(Data(spikes))/((epoch_times(epoch,2)- epoch_times(epoch,1))/10000);
            imagesc(NTC');h = colorbar;colormap(1-gray)
            set(h,'FontSize',text_font_size)
            title([ tfile_name ' PF Period: ' num2str(epoch) ' Rate: ' num2str(fq) 'Hz'],'FontSize',text_font_size)
            axis ij
            axis off
            subplot(n_epochs*2,2,subplot_count)
            subplot_count = subplot_count + 1;
            [SFx, SFy] = ScatterFields({spikes}, x_pos{epoch},y_pos{epoch});
            x = Data(x_pos{epoch});
            y = Data(y_pos{epoch});
            plot(x(1:10:end),y(1:10:end),'c.')
            title([ tfile_name ' ScatterFields Period: ' num2str(epoch)],'FontSize',text_font_size)
            axis ij
            axis tight
            axis off
            hold on
            plot(Data(SFx),Data(SFy),'rx')
            drawnow

            subplot(n_epochs*2,2,subplot_count:subplot_count+1)
            subplot_count = subplot_count + 2;
            r = Range(x_pos{epoch},'sec');
            plot(r(1:10:end),x(1:10:end)+y(1:10:end),'c')
            hold on
            %plot(r(1:10:end),y(1:10:end),'g.')
            plot(Range(SFx,'sec'),Data(SFx)+Data(SFy),'rx')
%            plot(Range(SFy,'sec'),Data(SFy),'kx')
            xlabel('sec','FontSize',text_font_size)
            ylabel('position (x + y)','FontSize',text_font_size)
            title([ tfile_name ' Position Vs. Time Period: ' num2str(epoch) ' Rate: ' num2str(fq) 'Hz'],'FontSize',text_font_size)
            %title('Position over time','FontSize',text_font_size)
            axis tight
            set(gca,'FontSize',text_font_size)
        end
        if save_figures
            saveas(f(figure_count-1),[tfile_name '_PF'],fig_string);
            close(f(figure_count-1))
        end
        
    end
    
end


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
