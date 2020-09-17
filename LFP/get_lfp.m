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
% ryan h 2020

p = inputParser;
p.addParameter('Fold',32000); % raw sample rate
p.addParameter('Fnew',1250); % downsampled sample rate for .lfp file
p.addParameter('overwrite_lfp',0); % reload wide band and overwrite .lfp 
p.addParameter('fs_for_datastruct',200); % downsampled sample rate for data struct
p.parse(varargin{:});
Fold = p.Results.Fold;
Fnew = p.Results.Fnew;
overwrite_lfp = p.Results.overwrite_lfp;
fs_for_datastruct = p.Results.fs_for_datastruct;

% if no inputs, we can still save the xml and .lfp files 
if ~exist('data','var')
    data.session_path = pwd;
    [~,data.basename] = fileparts(data.session_path);
end

warning off
% load lfp from .lfp file
if exist(fullfile(data.session_path,[data.basename,'.lfp']),'file') && overwrite_lfp == 0
    data = load_lfp_from_file(data,overwrite_lfp,fs_for_datastruct,Fnew,Fold);
    return
end

lfpfile = locate_ncs(data.session_path);

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

disp('loading downsampled lfp')
data = load_lfp_from_file(data,overwrite_lfp,fs_for_datastruct,Fnew,Fold);

end

function lfpfile = locate_ncs(session_path)
channels=table2cell(struct2table(dir([session_path,filesep,'*.ncs'])));
lfpfile=strcat(channels(:,2),filesep,channels(:,1));
[~,idx]=sort(str2double(extractBetween(lfpfile,'CSC','.ncs')));
lfpfile = lfpfile(idx);
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
if ~exist(probe_map{1},'file')
    msg = [probe_map{1},' does not exist'];
    error(msg)
end
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
    bz_MakeXMLFromProbeMaps(probe_map,data.session_path,basename,1,defaults)
end
session_info = LoadParameters(data.session_path);
end

function data = load_lfp_from_file(data,overwrite_lfp,fs_for_datastruct,Fnew,Fold)
    % get channel info
    lfpfile = locate_ncs(data.session_path);
    [probe_map,csc_list]=get_channel_list(lfpfile,data);
    session_info = make_load_xml(data,probe_map,overwrite_lfp,Fnew);
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