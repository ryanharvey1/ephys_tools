% PlotSpikes_script

anaDir=pwd;     % directory in which analyzed data is stored

load runpos startmaze endmaze X_run Y_run Vhs_cm startrun stoprun tracklength

cmbinsize=2;  % cm into which to bin the track
rereadcells=1;

cd(anaDir)
if ~rereadcells && exist('S_run.mat')
    load S_run S_run cellID
    nCells=length(S_run);
else
    tfiles = FindFiles('*.t');  % looks in folder you are in and all subfolders too
    S = LoadSpikes(tfiles);
    nCells=length(S);
    
    cellID=cell(nCells, 2);  % cell array with the ID of each cell, including tetrode number, and cell number on that tetrode
    S_maze=cell(nCells, 1);
    S_run=cell(nCells, 1);
    for c=1:nCells
        % restrict S to just maze epoch
        S_maze{c} = Restrict(S{c}, startmazes, endmazes);  % Restrict cell spikes to lap epoch time range
        S_run{c} = Restrict(S_maze{c}, startrun, stoprun);
        % find tetrode and cell number and store in cellID array
        [pathstr,name,ext]=fileparts(tfiles{c});
        cellID{c}(1)=str2num(name(3:4));
        perid=findstr(name, '.');
        if isempty(perid)
            perid=findstr(name, '_');
        end
        cellID{c}(2)=str2num(name(perid(end)+1:end));
    end
    save S_run S S_maze S_run cellID
    clear S S_maze
end

% load EEG from theta reference (should be the only CSC file in directory),
% filter for theta and find theta phase of all spikes
if exist('theta.mat')
    load theta eeg_theta
else
    eegfiles=findfiles('CSC*');
    [d, fname, ext] = fileparts(eegfiles{1});
    eeg_tsd = ReadCR_tsd([fname ext]);
    eeg_theta = Filter4Theta(eeg_tsd,6,11);  clear eeg_tsd;
    save theta eeg_theta;
end
[PH,thetapeaks] = ThetaPhase(S_run,eeg_theta,startmazes,endmazes);

rtracklength=round(tracklength/cmbinsize)*cmbinsize;
bins=-rtracklength:cmbinsize:rtracklength;
nBins=length(bins);
% find spike count in position bins, and occupancy
[FRonHStrack,Occ]=TuningCurvesBound(S_run, Vhs_cm, nBins, -rtracklength, rtracklength);

uisave({'FRonHStrack', 'Occ', 'cmbinsize', 'rtracklength'}, 'FRonHStrack3.mat')

% plot the positions in which spikes occurred (for each cell)
x_spikes=cell(nCells,1);
y_spikes=cell(nCells,1);
Vhs_spikes=cell(nCells,1);
for c=1:nCells
    % close figures so that matlab doesn't crash
    if exist('fh1', 'var')
        close(fh1)
    end
    if exist('fh2', 'var')
        close(fh2)
    end
    % plot the rat's path
    fh1=figure; plot(Data(X_run), Data(Y_run), 'k')
    title(['TT ', num2str(cellID{c}(1)), ' c ', num2str(cellID{c}(2))])
    xlabel('X (cm)')
    ylabel('Y (cm)')
    axis equal
    
    % interpolate to find the postions of each spike
    x_spikes{c} = interp1(Range(X_run), Data(X_run), Range(S_run{c}), 'linear');
    y_spikes{c} = interp1(Range(Y_run), Data(Y_run), Range(S_run{c}), 'linear');
    Vhs_spikes{c} = interp1(Range(Vhs_cm), Data(Vhs_cm), Range(S_run{c}), 'linear');
    % could also use [x_spikes, y_spikes, Vhs_spikes]=ScatterFields(S_maze, X_maze, Y_maze, Vhs)
    
    % plot each spike individually
    hold on; plot(x_spikes{c}, y_spikes{c}, 'r.')
    
    % main fig window for position plot (sized based on screen size)
    sc = [0 0 1 1];
    fh2 = figure;
    set(fh2,'Units','normalized','Position',[sc(1)+.05,sc(2)+.1,sc(3)*0.9,sc(4)*0.8]);
    orient tall;
    
    % make occupancy plot
    subplot(6,1,1); plot(bins, Occ)
    axis tight
    ylabel('Occupancy')
    title(['TT ', num2str(cellID{c}(1)), ' c ', num2str(cellID{c}(2))])
    
    % firing rate in postion bins plot
    subplot(6,1,2); plot(bins, FRonHStrack{c}./Occ)
    axis tight
    ylabel('Normalized firing rate')
    
    % make spikes by lap plot
    subplot(3,1,2); plot(Data(Vhs_cm), Range(Vhs_cm,'sec')/60, 'color', [0.6 0.6 0.6])
    hold on; plot(Vhs_spikes{c}, Range(S_run{c})/10000/60, 'r.')
    ylabel('Time on maze (min)')
    axis([-tracklength tracklength startmazes(1)/10000/60 endmazes(end)/10000/60])
    axis ij
    
    % make theta phase precession plot
    phase=Data(PH{c});
    if ~isempty(phase)
        subplot(3,1,3); plot(Vhs_spikes{c},phase(:,1),'k.', 'markersize', 2)
        hold on; plot(Vhs_spikes{c},phase(:,1)+1,'b.', 'markersize', 2)
    end
    xlabel('Position on track (deg)')
    ylabel('Theta phase')
    axis([-tracklength tracklength 0 2])
    axis xy
    
    % save the figures
    fnstrg = ['MazePlotR_' 'TT', num2str(cellID{c}(1)), 'c', num2str(cellID{c}(2))]; 
    saveas(fh1,[fnstrg '.tif']);
    
    fnstrg = ['ExpandedHorseshoePlotR_' 'TT', num2str(cellID{c}(1)), 'c', num2str(cellID{c}(2)), 'b']; 
    saveas(fh2,[fnstrg '.tif']);
end

