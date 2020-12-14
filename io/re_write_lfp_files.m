% re-write all .lfp files to 1250hz

% find all lfp files
% lfp_files = dir('F:\ClarkP30_Recordings\Data\**\*.lfp');
lfp_files = dir('F:\ClarkP30_Recordings\Data\**\*.nvt'); % loop through all folders

Fold = 32000;
Fnew = 1250;

WaitMessage = parfor_wait(length(lfp_files),'Waitbar',true);
% loop through each file
for i = 1:length(lfp_files)
    [~,basename] = fileparts(lfp_files(i).folder);
    
    % load raw data
    lfpfile = locate_ncs(lfp_files(i).folder);
    
    % locate xml 
    if isempty(dir([lfp_files(i).folder,filesep,'**\*.xml']))
        probe_map = get_channel_list(lfpfile,lfp_files(i).folder);
        write_xml(lfp_files(i).folder,probe_map,Fold,Fnew)
    end
    
    % preallocate & get time stamps
    [ts] = Nlx2MatCSC(lfpfile{1},[1 0 0 0 0], 0, 1, [] );
    signal = zeros(length(lfpfile),(length(ts)*512)*Fnew/Fold);
    
    % resample
    [fold, fnew] = rat(Fold./Fnew);
    for ii=1:length(lfpfile)
        fprintf(' %d ', ii);
        
        [~, filename] = fileparts(lfpfile{ii});
        ch = str2double(extractAfter(filename,'CSC'));
        
        [Samples]= Nlx2MatCSC(lfpfile{ii}, [0 0 0 0 1], 0, 1);
        
        signal(ch,1:length(Samples(:))*fnew/fold) =...
            resample(Samples(:), fnew, fold);
    end
    
    % save .lfp
    disp('saving .lfp file');
    fidout = fopen(fullfile(lfp_files(i).folder,[basename,'.lfp']), 'w');
    fwrite(fidout,signal,'int16');
    fclose(fidout);
    WaitMessage.Send;
    
end
WaitMessage.Destroy;

function lfpfile = locate_ncs(session_path)
channels=table2cell(struct2table(dir([session_path,filesep,'*.ncs'])));
lfpfile=strcat(channels(:,2),filesep,channels(:,1));
[~,idx]=sort(str2double(extractBetween(lfpfile,'CSC','.ncs')));
lfpfile = lfpfile(idx);
end

function write_xml(path,map,Fold,FNew)
% Wrapper for buzcode bz_MakeXMLFromProbeMaps with defaults

defaults.NumberOfChannels = 1;
defaults.SampleRate = Fold;
defaults.BitsPerSample = 16;
defaults.VoltageRange = 20;
defaults.Amplification = 1000;
defaults.LfpSampleRate = FNew;
defaults.PointsPerWaveform = 64;
defaults.PeakPointInWaveform = 32;
defaults.FeaturesPerWave = 4;
[~,basename] = fileparts(path);
bz_MakeXMLFromProbeMaps({map},path,basename,1,defaults)
end

function probe_map = get_channel_list(lfpfile,path)

% many session may only have one continuous 'lfp' channel per tetrode
if length(dir(fullfile(path,'*.ncs'))) ==...
        length(dir(fullfile(path,'*.ntt')))
    probe_map = {[num2str(max(str2double(extractBetween(lfpfile,'CSC','.ncs')))),...
        'tt_without_all_channels.xlsx']};
else
    probe_map = {[num2str(max(str2double(extractBetween(lfpfile,'CSC','.ncs')))/4),...
        'tt.xlsx']};
end
end