function filter_raw_dat(varargin)
% filter_raw_dat: high pass filters raw neuralynx CSC files and writes them
% to 'filtered.dat'
%
% Ryan Harvey 2019

% parse input
p = inputParser;
p.addParameter('output_file','filtered.dat');
p.addParameter('raw_dir',pwd); % folder with CSCs
p.addParameter('hi_pass',300);
p.addParameter('fs',32000);
p.addParameter('verbose',true);

p.parse(varargin{:});

output_file = ['temp_',p.Results.output_file];
raw_dir = p.Results.raw_dir;
hi_pass = p.Results.hi_pass;
fs = p.Results.fs;
verbose = p.Results.verbose;

cd(raw_dir)

if verbose
    tic;
    disp('Making hi-pass filtered dat file from CSCs')
end

%locate lfp files
chs=dir('*.ncs');

% open filtered dat file for writing
fid_filt = fopen(output_file, 'w');

% create high pass filter
[b1, a1] = butter(3, hi_pass/fs*2, 'high');

% loop though channels
for ch=1:length(chs)
    
    if verbose
        disp(['CSC',num2str(ch),'.ncs'])
    end
    
    % extract CSC data with mex
    if ismac==1
        [dataRAW]= Nlx2MatCSC_v3(['CSC',num2str(ch),'.ncs'], [0 0 0 0 1], 0, 1);
    else
        [dataRAW]= Nlx2MatCSC(['CSC',num2str(ch),'.ncs'], [0 0 0 0 1], 0, 1);
    end
    
    % high pass filter data with forwards backwards filter
    datr = filtfilt(b1, a1, dataRAW(:));
    
    % gather data and convert to int16
    datcpu  = gather_try(int16(datr));
    
    % write to disk
    fwrite(fid_filt, datcpu, 'int16');
end
% close file
fclose(fid_filt);

clear datcpu dataRAW datr 




% This final section transpose binary file (this results in a 300 %
% increase in speed when pulling out waveforms later on)
filenamestruct = dir(output_file);
dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(ch*dataTypeNBytes);  % Number of samples per channel
mmf = memmapfile(output_file, 'Format', {'int16', [nSamp ch], 'x'});
x=mmf.Data.x;

fid_filt = fopen(erase(output_file,'temp_'), 'w');

% create batches
batch=ceil(linspace(0,nSamp,ceil(nSamp/fs/4)));

% loop though batches
for i=1:length(batch)-1
    
    if verbose
        disp(['batch ',num2str(batch(i)+1),' to ',num2str(batch(i+1)),...
            '   ',num2str(i),' of ',num2str(length(batch)-1)])
    end
    
    % pull out batch of data
    datr = x(batch(i)+1:batch(i+1),1:ch);
    
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
end