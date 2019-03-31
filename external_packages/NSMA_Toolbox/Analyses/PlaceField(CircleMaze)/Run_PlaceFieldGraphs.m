function Run_PlaceFieldGraphs(event_ts1, event_ts2, horseshoe_flag)

% Run_PlaceFieldGraphs  Produces lap-by-lap firing rasters and phase precession plots
% 
% Run_PlaceFieldGraphs(event_ts1, event_ts2, horseshoe_flag)
% 
% INPUTS:
%   event_ts1 & 2 - array of timestamps of events for up to 2 different types of events.  These
%       are optional - may omit these inputs, or just input one type of event (_ts1) but not _ts2).
%       If horseshoe_flag input is specified, need to enter something for both events (may be empty lists []).
%   horseshoe_flag - 1=animal runs in "horseshoe" pattern on circular track, where it runs once
%       around circle in one direction, turns and runs once around in reverse direction, etc.
%       0(or no inputs)=animal runs regular laps around circle
% OUTPUTS:
%   (none)
%
% User runs program by being shown a lap-by-lap firing raster ("SpotPlot") for each cell, and 
%   asked to click in center of place field.  Grey vertical lines represent spikes in each lap,
%   and circles represent firing rate - bigger diameter=higher rate.  After clicking on center
%   of place field, user is asked whether there is another place field for that cell.  If yes,
%   user types "y" [enter], at which point (s)he may click on another field.  If no, just hit
%   [enter], and the graph for the next cell pops up.
% After this interactive portion of the program, 3 types of plots are produced for each cell and
%   Matlab .fig files are saved in current folder:
%   1) "SpotPlot" - lap-by-lap raster showing positions of spikes and possibly other events.
%       Larger diameter of circles = higher firing rate.
%   2) "PhaseTime" - lap-by-lap raster with "theta time" (theta cycles prior to and after rat 
%       runs through center of place field) rather than position on x-axis
%   3) Phase precession plot - shows phase precession of spikes as rat passes through center of
%       place field
%
% PL, modified by MN spr.'02

cd F:\Users\BClarkLab\Desktop\Ryan\Test_session\2016-10-27_18-10-02p
if nargin < 3
    horseshoe_flag = 0;
end
if nargin < 2
    event_ts2 = [];
end
if nargin < 1
    event_ts1 = [];
end


% get position data
if exist('position.mat')
    load position X Y;
    ts_range = Range(X,'ts');
    ts_MazeStart = ts_range(1);
    ts_MazeEnd = ts_range(end);
    [V,W] = ParameterizeToCircularTrack(X,Y);
end %if exist('position.mat')


% get spike data
TTfolders = dir;
tfiles = {};
S = {};
for i=3:length(TTfolders)
    if TTfolders(i).isdir
        if ( strcmp(TTfolders(i).name(1:2),'TT') | strcmp(TTfolders(i).name(1:2),'tt') )
            videofile=textread(cell2mat(FindFiles('VT1.txt'))); % Ryan H 12/9/16
            cd(TTfolders(i).name)
            
            tfiles_temp = FindFiles('*.t'); 
            S_temp = LoadSpikes(tfiles_temp);

            tfiles = [tfiles; tfiles_temp];
            S = [S; S_temp];
            
            cd ..
            
        end %if ( strcmp ...
    end %if TTfolders ...
end %for
if ~exist('position.mat') % Ryan Harvey 12/9/16
    ts_MazeStart = videofile(1);
    ts_MazeEnd = videofile(end,1); 
    V=videofile(:,2);
end
% Get eeg data
if exist('theta.mat')
    load theta eeg_theta
else
    eegfiles=FindFiles('CSC*');
    %eegfiles=findfiles('*.Ncs');
    [d, fname, ext] = fileparts(eegfiles{1});
    eeg_tsd = ReadCR_tsd([fname ext]);
    eeg_theta = Filter4Theta(eeg_tsd,6,11);  clear eeg_tsd;
end
PH = ThetaPhase(S,eeg_theta,ts_MazeStart,ts_MazeEnd);


% Calculate couple other variables used in plot functions
AmplifierGainFactor = 2000;
A2DGainFactor = 2;
eeg_scalefactor = 10^6 * 4 /(4096*AmplifierGainFactor*A2DGainFactor);

% Vts = Range(V,'ts'); Ryan Harvey 12/9/16
 Vts = range(V,'includenan');
nCells = length(S);


% Interactive part where user clicks on place fields
graph = 1;
pf_center = []; cellnum = []; field = [];
% loop over cells
for iC = 1:nCells
    
    cell_ts = Data(S{iC});
    
    again = 1;
    iF = 0;
    while again
        iF = iF+1;

        [Ddir, fname, ext] = fileparts(tfiles{iC});
        
        % convert underscores to dashes and / to //
        fname = strrep(fname,'_','-');
        fname = strrep(fname,'/','//');
        
        SpotPlot_quick(V, cell_ts, fname, horseshoe_flag, event_ts1, event_ts2);
        
        disp('Click in center of place field, or if no field, click anywhere on plot.');
    
        [pfcenter, ydum] = ginput(1);
        pf_center = [pf_center; pfcenter];
        cellnum = [cellnum; iC];
        field = [field; iF];
        graph = graph+1;

        rstr = input('Another place field on same cell (y(es)/[n](o))?: ','s');
        if strcmpi(rstr,'y')
            again = 1;
            close all;
        else
            again = 0;
        end
        
    end% while
    
    close all;
      
end %for iC

save pf_vals cellnum field pf_center;
close all;


% Produce place field graphs
for igraph = 1:length(cellnum) 
    
    iC = cellnum(igraph);
    iF = field(igraph);
    pfc = pf_center(igraph);
    
    cell_ts = Data(S{iC});

    [Ddir, fname, ext] = fileparts(tfiles{iC});

    titlestr = [fname ' - Field ' num2str(iF)]; 
    pow = Data(eeg_theta).^2;
    ts_eeg  = Range(eeg_theta,'ts');
    [ts_env,pow_env, PeaksIdx] = PowerEnvelope(ts_eeg,pow);
    thetapow_tsd = tsd(ts_env,pow_env);
    fh1 = PlaceFieldSpotPlot(V,cell_ts,horseshoe_flag,event_ts1,event_ts2,pfc,titlestr,thetapow_tsd);
    fh2 = PhaseTimePlot(V,PH{iC},eeg_theta,horseshoe_flag,event_ts1,event_ts2,pfc,15,titlestr);
    fh3 = PhasePrecessPlot(V,PH{iC},eeg_theta,horseshoe_flag,pfc,50,eeg_scalefactor,titlestr);

    
    % Save graphs
    fnstrg = ['SpotPlot_' fname '_Field' num2str(iF)]; 
    saveas(fh1,[fnstrg '.fig']);

    fnstrg = ['Phasetime_' fname '_Field' num2str(iF)]; 
    saveas(fh2,[fnstrg '.fig']);

    fnstrg = ['PhasePrecess_' fname '_Field' num2str(iF)]; 
    saveas(fh3,[fnstrg '.fig']);
    
    
    close all;

end % for igraph