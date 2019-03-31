% Kilosort_script.m
%
% Runs kilosort on multiple days, writes .spikes and .clu files for TTs and
% writes spike times (.t files) for Si probes.
%
% Data directories should look like this:
%
%   pathIn\
%       R###\
%           R###-YYYY-MM-DD\
%               raw-intan-files.dat
%               ...
%
%   pathOut\
%       R###\
%           R###_channel_map_24TT.txt
%           R###_channel_map_32Si.txt
%           R###-YYYY-MM-DD\
%               sd file, video tracking, keys files, etc...
%               will write .t files here
%
% Aug 2017
% Brendan Hasz
% haszx010@umn.edu
% David Redish Lab, University of Minnesota Twin Cities

%% Settings to run Kilosort

% Input and output data directories
pathIn = 'F:\Data\PrePromoted';
pathOut = 'F:\Data\Promoted';

% Sensor to use (24TT or 32Si)
sensors = {'24TT', '32Si'}; %do both TTs and probes for all days
%sensors = {'24TT'};
%sensors = {'32Si'};

% Configuration file prefix (will append sensor.m, e.g.: "config_KiloSort_24TT.m")
config_prefix = 'C:\Users\hasz\Documents\MATLAB\BMHcodeset\Neural\KiloSort\IntanToKilosort\config_KiloSort_';

% List of Rat/Days to process
ssns = {'R422-2017-06-07'};
%{
ssns = {...
    'R422-2017-05-30', ...
    'R422-2017-05-31', ...
    'R422-2017-06-01', ...
    'R422-2017-06-02', ...
    'R422-2017-06-03', ...
    'R422-2017-06-04', ...
    'R422-2017-06-05', ...
    'R422-2017-06-06', ...
    'R422-2017-06-07' ...
    };
%}


%% Run kilosort on each day and sensor
for iD = 1:length(ssns)
    ssn = ssns{iD};
    for iS = 1:length(sensors)
        sensor = sensors{iS};
        
        fprintf('\n\nRunning KiloSort on %s %s...\n', ssn, sensor)
        
        % Get file/path names
        rat = ssn(1:4);
        fpathIn = [pathIn filesep rat filesep ssn];
        fpath = [pathOut filesep rat filesep ssn filesep 'KiloSort-' sensor];
        data_fname = [ssn '-KilosortRaw-' sensor '.dat'];
        chanMap_fname = [pathOut filesep rat filesep rat '_channel_map_' sensor '.txt'];
        config_fname = [config_prefix sensor '.m'];
        
        % Make output directory if it doesn't exist
        if ~exist(fpath, 'dir')
            mkdir(fpath);
        end

        % Start log file
        fname = sprintf('%s-%s-Kilosort-%s.log', ...
            ssn, sensor, datestr(now, 'yyyy-mm-dd-HH-MM-SS'));
        fid = fopen(fname,'w');
        fprintf(fid,'%s %s\n', ssn, sensor);
        fprintf(fid,'Input path: %s\n', pathIn);
        fprintf(fid,'Output path: %s\n', pathOut);

        % Run Kilosort
        fprintf(fid,'Kilosort started at %s\n',datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        tic;
        run_KiloSort(fpathIn, fpath, data_fname, chanMap_fname, config_fname);
        ttr = toc; toc

        % Close log file
        fprintf(fid,'Kilosort finished at %s\n',datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        fprintf(fid,'Elapsed time: %s\n', num2str(ttr));
        fclose(fid);
    end
end


%% Now sort in Phy

% Manually.  Wheeee!
warning('SORT CLUSTERS IN PHY BEFORE PROCEEDING!')
pause;


%% Write .spikes and .clu files for TTs

sensor = '24TT';
ssns = {'R395-2017-03-13', 'R395-2017-03-14'};

for iD = 1:length(ssns)
    ssn = ssns{iD};
    rat = ssn(1:4);
    timePath = [pathIn filesep rat filesep ssn];
    fpath = [pathOut filesep rat filesep ssn];
    fpathIn = [fpath filesep 'KiloSort-' sensor];
    rawFname = [fpathIn filesep ssn '-KilosortRaw-' sensor '.dat'];
    chanMapFname = [pathOut filesep rat filesep rat '_channel_map_' sensor '.txt'];

    % Extract the spikes and save .spikes and .clu files
    % (this could take 2-3 hrs per session)
    ExtractSpikesFromKilosort(rawFname, timePath, fpathIn, fpath, chanMapFname);
end


%% Write spike times for probes

ssns = {'R395-2017-04-02'};
sensor = '32Si';

for iD = 1:length(ssns)
    ssn = ssns{iD};
    disp(ssn)
    rat = ssn(1:4);
    timePath = [pathIn filesep rat filesep ssn];
    fpath = [pathOut filesep rat filesep ssn];
    fpathIn = [fpath filesep 'KiloSort-' sensor];
    SaveSpikeTimesFromKilosort(timePath, fpathIn, fpath, sensor);
end