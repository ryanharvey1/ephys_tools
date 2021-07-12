function filter_raw_dat_from_dat(varargin)
% filter_raw_dat: high pass filters raw neuralynx CSC files and writes them
% to 'filtered.dat'
%
% Ryan Harvey 2019
% Modified for probes Laura Berkowitz 2020
addpath(genpath('D:\Users\BClarkLab\ephys_tools\external_packages\analysis-tools'))

% parse input
p = inputParser;
p.addParameter('output_file','filtered.dat');
p.addParameter('raw_dir',pwd); % folder with CSCs
p.addParameter('hi_pass',300);
p.addParameter('fs',30000);
p.addParameter('verbose',true);
p.addParameter('nChannels',64);

p.parse(varargin{:});

output_file = ['temp_',p.Results.output_file];
raw_dir = p.Results.raw_dir;
hi_pass = p.Results.hi_pass;
fs = p.Results.fs;
verbose = p.Results.verbose;
nChannels = p.Results.nChannels;

session_path = raw_dir;
[~,basename] = fileparts(session_path);

if verbose
    tic;
    disp('Making hi-pass filtered dat file from memmory mapped dat')
end

if ~exist(fullfile(session_path,[basename,'.dat']),'file')
    % Get json path
    jsonfile = locate_binary(session_path);
    
    % Load mapped data
    D = load_open_ephys_binary(jsonfile{1},'continuous',1,'mmap');
    
    % Save headstage data and accelerometer data in separate dat files
    data_idx = contains({D.Header.channels.description},'Auxiliar');
    
    % Open filtered dat file for writing
    fid_filt = fopen(output_file, 'w');
    
    % Create high pass filter
    [b1, a1] = butter(3, hi_pass/fs*2, 'high');
    
    % Filter each channel
    fprintf('filtering channel...');
    ch = find(~data_idx');
    for ii = ch'
        fprintf(' %d ', ii);
        
        % high pass filter data with forwards backwards filter
        datr = filtfilt(b1, a1, double(D.Data.data.mapped(ii,:)));
        
        % gather data and convert to int16
        datcpu  = gather_try(int16(datr));
        
        % write to disk
        fwrite(fid_filt, datcpu, 'int16');
    end
    
    % close file
    fclose(fid_filt);
    
    clear datcpu dataRAW datr
else
    % Load .dat
    fileName = fullfile(session_path,[basename,'.dat']);
    filenamestruct = dir(fileName);
    dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8')); % determine number of bytes per sample
    nSamp = filenamestruct.bytes/(nChannels*dataTypeNBytes);  % Number of samples per channel
    mmf = memmapfile(fileName, 'Format', {'int16', [nChannels nSamp], 'Data'});
    % Open filtered dat file for writing
    fid_filt = fopen(output_file, 'w');
    
    % Create high pass filter
    [b1, a1] = butter(3, hi_pass/fs*2, 'high');
    
    % Filter each channel
    fprintf('filtering channel...');
    for ii = 1:nChannels
        fprintf(' %d ', ii);
        
        % high pass filter data with forwards backwards filter
        datr = filtfilt(b1, a1, double(mmf.Data.Data(ii,:)));
        
        % gather data and convert to int16
        datcpu  = gather_try(int16(datr));
        
        % write to disk
        fwrite(fid_filt, datcpu, 'int16');
    end
    
    % close file
    fclose(fid_filt);
    
    clear datcpu dataRAW datr
end


% This final section transpose binary file (this results in a 3x
% increase in speed when pulling out waveforms later on)
filenamestruct = dir(output_file);
dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(nChannels*dataTypeNBytes);  % Number of samples per channel
mmf = memmapfile(output_file, 'Format', {'int16', [nSamp nChannels], 'x'});
x=mmf.Data.x;

fid_filt = fopen(erase(output_file,'temp_'), 'w');

% create batches
batch = ceil(linspace(0,nSamp,ceil(nSamp/fs/4)));

% loop though batches
for i=1:length(batch)-1

    if verbose
        disp(['batch ',num2str(batch(i)+1),' to ',num2str(batch(i+1)),...
            '   ',num2str(i),' of ',num2str(length(batch)-1)])
    end

    % pull out batch of data
    datr = x(batch(i)+1:batch(i+1),1:nChannels);

    % gather data and convert to int16 & transpose
    datcpu = gather_try(int16(datr'));

    % write to disk
    fwrite(fid_filt, datcpu, 'int16');
end
% close file
fclose(fid_filt);

disp('Cleaning up temp file')

clear x mmf

delete(output_file)

if verbose
    toc
end

% Local Functions below
    function lfpfile = locate_binary(session_path)
        file = table2cell(struct2table(dir([session_path,filesep,'**\*.oebin'])));
        lfpfile = strcat(file(:,2),filesep,file(:,1));
    end

end