function data = get_lfp(data,varargin)
% get_lfp: loads Nlx csc lfp data, resamples it, and gets a few extra theta features
%
% resample applies an antialiasing FIR lowpass filter to the signal and
% compensates for the delay introduced by the filter
%
%   Input:
%           data: ephys_tool's data structure
%           varargin:
%               Fold: lfp sample rate (32000)
%               Fnew: resampled lfp sample rate (1000)
%               overwrite_lfp: overwrite the .lfp and xml (0)
%               fs_for_datastruct: resampled lfp sample rate to be saved in
%               the postprocessed data structure (200)
%
% You can also call get_lfp without the 'data' input. Just make sure to cd to your
% raw data folder beforehand.
%
%   Output:
%           data:
%               data.lfp.theta_phase
%               data.lfp.ts
%               data.lfp.signal
%               data.lfp.lfpsamplerate
%               data.lfp.load_date
%               data.lfp.creation_date
%               data.lfp.csc_channel_list
%
% ryan h 2020, adapted for probes laura b Nov 2020
addpath(genpath('D:\Users\BClarkLab\ephys_tools\external_packages\analysis-tools'))
p = inputParser;
p.addParameter('Fold',32000); % raw sample rate
p.addParameter('Fnew',1250); % downsampled sample rate for .lfp file
p.addParameter('overwrite_lfp',0); % reload wide band and overwrite .lfp
p.addParameter('fs_for_datastruct',200); % downsampled sample rate for data struct
p.addParameter('probe_map','channel_map_64.xlsx') % default for data already mapped to disk

p.parse(varargin{:});
Fold = p.Results.Fold;
Fnew = p.Results.Fnew;
overwrite_lfp = p.Results.overwrite_lfp;
fs_for_datastruct = p.Results.fs_for_datastruct;
probe_map = p.Results.probe_map;


% if no inputs, we can still save the xml and .lfp files
if ~exist('data','var')
    data.session_path = pwd;
    [~,data.basename] = fileparts(data.session_path);
end

% Check for nsc files
lfpfile = locate_ncs(data.session_path);

warning off
% load lfp from .lfp file
if exist(fullfile(data.session_path,[data.basename,'.lfp']),'file') && overwrite_lfp == 0 && ~isempty(lfpfile)
    data = load_lfp_from_file(data,overwrite_lfp,fs_for_datastruct,Fnew,Fold);
    return
elseif exist(fullfile(data.session_path,[data.basename,'.lfp']),'file') && overwrite_lfp == 0 && isempty(lfpfile)
    Fold = 30000;
    data = load_lfp_created_from_dat(data,overwrite_lfp,fs_for_datastruct,Fnew,Fold,probe_map);
    return
end



% Load lfp from Neuralynx ncs files
if ~isempty(lfpfile)
    
    % preallocate & get time stamps
    [ts] = Nlx2MatCSC(lfpfile{1},[1 0 0 0 0], 0, 1, [] );
    signal = zeros(length(lfpfile),(length(ts)*512)*Fnew/Fold);
    
    %account for non-integer fs
    % get new fs exact
    [fold, fnew] = rat(Fold./Fnew);
    
    % loop though each channel and resample lfp
    fprintf('channel...');
    for ii=1:length(lfpfile)
        fprintf(' %d ', ii);
        
        [~, filename] = fileparts(lfpfile{ii});
        ch = str2double(extractAfter(filename,'CSC'));
        
        [Samples]= Nlx2MatCSC(lfpfile{ii}, [0 0 0 0 1], 0, 1);
        
        signal(ch,1:length(Samples(:))*fnew/fold) = resample(Samples(:), fnew, fold);
    end
    fprintf('lfp loaded\n');
    
    % save to folder for quick reading
    disp('saving .lfp file');
    fidout = fopen(fullfile(data.session_path,[data.basename ,'.lfp']), 'w');
    fwrite(fidout,signal,'int16');
    fclose(fidout);
    
    % Load 200Hz LFP
    disp('loading downsampled lfp')
    data = load_lfp_from_file(data,overwrite_lfp,fs_for_datastruct,Fnew,Fold);
    
    % Load lfp from dat file
else
    % Set some basic 
    Fold = 30000; % Open ephys sampling rate
    nChannels = length(xlsread(probe_map));
    
    % Get json path
    jsonfile = locate_binary(data.session_path);
    
    % Load mapped data
    D = load_open_ephys_binary(jsonfile{1},'continuous',1,'mmap');
    
    if ~exist(fullfile(data.session_path,[data.basename,'_accel.dat']),'file')
        % Save headstage data and accelerometer data in separate dat files
        data_idx = contains({D.Header.channels.description},'Auxiliar');
        disp('saving separate .dat files for headstage data and accelerometer data');
        
        % headstage data
        fidout = fopen(fullfile(data.session_path,[data.basename ,'.dat']), 'w');
        fwrite(fidout,D.Data.Data.mapped(~data_idx',:),'int16');
        fclose(fidout);
        
        % accelerometer data
        fidout = fopen(fullfile(data.session_path,[data.basename ,'_accel.dat']), 'w');
        fwrite(fidout,D.Data.Data.mapped(data_idx',:),'int16');
        fclose(fidout);
        
        % Account for non-integer fs
        % get new fs exact
        [fold, fnew] = rat(Fold./Fnew);
        
        % loop though each channel and resample lfp
        fprintf('channel...');
        ch = find(~data_idx');
        for ii = ch'
            fprintf(' %d ', ii);
            signal(ii,:) = resample(double(D.Data.data.mapped(ii,:)), fnew, fold);
        end
        
        % save to folder for quick reading
        disp('saving .lfp file');
        fidout = fopen(fullfile(data.session_path,[data.basename ,'.lfp']), 'w');
        fwrite(fidout,signal,'int16');
        fclose(fidout);
    else 
        % Load .dat 
        fileName = fullfile(data.session_path,[data.basename,'.dat']);
        filenamestruct = dir(fileName);
        dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8')); % determine number of bytes per sample
        nSamp = filenamestruct.bytes/(nChannels*dataTypeNBytes);  % Number of samples per channel
        mmf = memmapfile(fileName, 'Format', {'int16', [nChannels nSamp], 'Data'});

       % Account for non-integer fs
        % get new fs exact
        [fold, fnew] = rat(Fold./Fnew);
        
        % loop though each channel and resample lfp
        fprintf('channel...');
        for ii = 1:nChannels
            fprintf(' %d ', ii);
            signal(ii,:) = resample(double(mmf.Data.Data(ii,:)), fnew, fold);
        end
        
        % save to folder for quick reading
        disp('saving .lfp file');
        fidout = fopen(fullfile(data.session_path,[data.basename ,'.lfp']), 'w');
        fwrite(fidout,signal,'int16');
        fclose(fidout); 
    end
    
    
    disp('loading downsampled lfp')
    data = load_lfp_created_from_dat(data,overwrite_lfp,fs_for_datastruct,Fnew,Fold,probe_map);
    
end



end

function lfpfile = locate_ncs(session_path)
channels=table2cell(struct2table(dir([session_path,filesep,'*.ncs'])));
lfpfile=strcat(channels(:,2),filesep,channels(:,1));
[~,idx]=sort(str2double(extractBetween(lfpfile,'CSC','.ncs')));
lfpfile = lfpfile(idx);
end

function lfpfile = locate_binary(session_path)
file = table2cell(struct2table(dir([session_path,filesep,'**\*.oebin'])));
lfpfile = strcat(file(:,2),filesep,file(:,1));
end

function [probe_map,csc_list]=get_channel_list(lfpfile,data)

% many session may only have one continuous 'lfp' channel per tetrode
if length(dir(fullfile(data.session_path,'*.ncs'))) ==...
        length(dir(fullfile(data.session_path,'*.ntt')))
    probe_map = {[num2str(max(str2double(extractBetween(lfpfile,'CSC','.ncs')))),...
        'tt_without_all_channels.xlsx']};
    csc_list.channel_num(:,1) = 1:4:max(str2double(extractBetween(lfpfile,'CSC','.ncs'))) * 4 - 3;
    csc_list.tetrode_num(:,1) = 1:max(str2double(extractBetween(lfpfile,'CSC','.ncs')));
else
    probe_map = {[num2str(max(str2double(extractBetween(lfpfile,'CSC','.ncs')))/4),...
        'tt.xlsx']};
    csc_list.channel_num(:,1) = str2double(extractBetween(lfpfile,'CSC','.ncs'));
    csc_list.tetrode_num(:,1) = reshape(repmat(1:length(lfpfile)/4, 4, 1),length(lfpfile)/4*4,1);
end
end

function session_info = make_load_xml(data,probe_map,overwrite_lfp,Fnew,Fold)
% if ~exist(probe_map{1},'file')
%     msg = [probe_map{1},' does not exist'];
%     error(msg)
% end
if ~exist([fullfile(data.session_path,data.basename),'.xml'],'file') || overwrite_lfp
    % write xml
    defaults.NumberOfChannels = 1;
    defaults.SampleRate = Fold;
    defaults.BitsPerSample = 16;
    defaults.VoltageRange = 20;
    defaults.Amplification = 1000;
    defaults.LfpSampleRate = Fnew;
    defaults.PointsPerWaveform = 64;
    defaults.PeakPointInWaveform = 32;
    defaults.FeaturesPerWave = 4;
    [~,basename] = fileparts(data.session_path);
    bz_MakeXMLFromProbeMaps({probe_map},data.session_path,basename,1,defaults)
end
session_info = LoadParameters(data.session_path);
end

function data = load_lfp_from_file(data,overwrite_lfp,fs_for_datastruct,Fnew,Fold)
% get channel info
lfpfile = locate_ncs(data.session_path);
[probe_map,csc_list]=get_channel_list(lfpfile,data);
session_info = make_load_xml(data,probe_map,overwrite_lfp,Fnew,Fold);
info = dir(fullfile(data.session_path,[data.basename,'.lfp']));

lfp = bz_GetLFP('all','basepath',data.session_path,...
    'basename',data.basename,...
    'noPrompts',true,...
    'downsample',1);

% resample to 200 hz or other sampling rate
[fold, fnew] = rat(lfp.samplingRate./fs_for_datastruct);
lfp.data = resample(double(lfp.data), fnew, fold);
lfp.samplingRate = fs_for_datastruct;

% mark good channels
good = zeros(1,length(lfpfile));
% check spike group from xml file
good([session_info.spikeGroups.groups{:}] + 1) = 1;
% check if a channel is all zeros
csc_list.good_channels(sum(lfp.data) ~= 0 & good,1) = 1;

% get aligned ts
if ~isfield(data,'offset')
    [vts] = Nlx2MatVT(fullfile(data.session_path,'VT1.nvt'),[1,0,0,0,0,0],0,1);
    data.offset = vts(1);
end
ts = Nlx2MatCSC(lfpfile{1},[1 0 0 0 0], 0, 1, [] );
ts = interp1(linspace(1,length(lfp.data),length(ts)), ts, 1:length(lfp.data));
ts = (ts-data.offset) / 10^6;

% get phase
[theta_phase,~,~]=Phase([ts',BandpassFilter(double(lfp.data),lfp.samplingRate, [4 12])]);

data.lfp.theta_phase = theta_phase(:,2:end)';
data.lfp.ts = ts;
data.lfp.signal=double(lfp.data)';
data.lfp.lfpsamplerate=lfp.samplingRate;
data.lfp.load_date = date;
data.lfp.creation_date = info.date;
data.lfp.channel_list = csc_list;

warning on
end

function data = load_lfp_created_from_dat(data,overwrite_lfp,fs_for_datastruct,Fnew,Fold,probe_map)
% get channel info
file = table2cell(struct2table(dir([data.session_path,filesep,'**\*.lfp'])));

% Assume Cambridge Neurotech E1 64 channel probe for now (channels must
% be remapped in Open Ephys).
make_load_xml(data,probe_map,overwrite_lfp,Fnew,Fold);
info = dir(fullfile(data.session_path,[data.basename,'.lfp']));

% LoadXml using KS2Wrapper function
d   = dir('*.xml');
par = LoadXml(fullfile(data.session_path,d(1).name));

lfp = bz_GetLFP('all','basepath',data.session_path,...
    'basename',data.basename,...
    'noPrompts',true,...
    'downsample',1);

% resample to 200 hz or other sampling rate
[fold, fnew] = rat(lfp.samplingRate./fs_for_datastruct);
lfp.data = resample(double(lfp.data), fnew, fold);
lfp.samplingRate = fs_for_datastruct;
nChannels = length(lfp.channels);

% initialize csc_list struct
csc_list.good_channels = zeros(nChannels,1);
% check spike group from xml file (bad channels are skipped)
good = ~logical([par.AnatGrps(:).Skip]);
% check if a channel is all zeros
csc_list.good_channels(sum(lfp.data) ~= 0 & good,1) = 1;
% phy outputs channels in 0 to n channels.
csc_list.channel_num = [1:nChannels]';
tetrode_num = [repmat([1:4]',1,nChannels/4)]';
csc_list.tetrode_num = tetrode_num(:);

% get resampled timestamps
ts = interp1(linspace(1,length(lfp.data),length(lfp.timestamps)), lfp.timestamps, 1:length(lfp.data));

% get phase
[theta_phase,~,~] = Phase([ts',BandpassFilter(double(lfp.data),lfp.samplingRate, [4 12])]);

data.offset = 0;
% Save resampled lfp back to file
data.lfp.theta_phase = theta_phase(:,2:end)';
data.lfp.ts = ts;
data.lfp.signal=double(lfp.data)';
data.lfp.lfpsamplerate=lfp.samplingRate;
data.lfp.load_date = date;
data.lfp.creation_date = info.date;
data.lfp.channel_list = csc_list;

warning on

end

