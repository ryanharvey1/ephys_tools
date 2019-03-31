function CopyDat(fnameIn,fnameOut,varargin)

% USAGE
%     CopyDat(fnameIn,fnameOut,'optionNames',optionalues)
% copy a given number of seconds of data from one dat file into another one    
%     
% INPUT    
%     - fnameIn     : source file
%     - fnameOIut   : target file
% options:
%     - 'duration'  : duration in seconds
%     - 'start'     : start time in seconds
%     - 'channelIx' : subset of channels to be copied

% Adrien Peyrache, 2012

start = 0;
channelIx = [];
duration = [];

for i = 1:2:length(varargin),
  if ~isa(varargin{i},'char'),
    error(['Parameter ' num2str(i+3) ' is not a property (type ''help CopyDat'' for details).']);
  end
  switch(lower(varargin{i})),
    case 'start',
      start = varargin{i+1};
      if ~isa(start,'double')
        error('Incorrect value for property ''start''. Should be char array.');
      end
    case 'channelix',
      channelIx = varargin{i+1};
      if ~isa(channelIx,'double')
        error('Incorrect value for property ''chanellIx''. Should be char array.');
      end
    case 'duration',
      duration = varargin{i+1};
      if ~isa(duration,'double')
        error('Incorrect value for property ''duration''. Should be char array.');
      end  
  end
end

fxml = [fnameIn(1:end-4) '.xml'];
if ~exist(fnameIn,'file') && ~exist(fxml,'file')
  error('Dat file and/or Xml file does not exist')
end
sizeInBytes = 2; % change it one day...
chunk = 1e5; % depends on the system... could be bigger I guess

syst = LoadXml(fxml);

fInfo = dir(fnameIn);

nbChan = syst.nChannels;
if ~isempty(duration)
    nBytes = syst.SampleRate*duration*sizeInBytes*nChan;
    if fInfo.bytes<nBytes
        error('Duration longer than total duration of file')
    end
else
    nBytes = fInfo.bytes;
end
nbChunks = floor(nBytes/(nbChan*sizeInBytes*chunk));

fidO = fopen(fnameOut,'w');
fidI = fopen(fnameIn,'r');

start = floor(start*syst.SampleRate)*syst.nChannels*sizeInBytes;
status = fseek(fidI,start,'bof');
if status ~= 0,
		error('Could not start reading (possible reasons include trying to read a closed file or past the end of the file).');
end

for ii=1:nbChunks
    h=waitbar(ii/(nbChunks+1));
    dat = fread(fidI,nbChan*chunk,'int16');
    if ~isempty(channelIx)
        dat = reshape(dat,[nbChan chunk]);
        dat = dat(channelIx,:);
        dat = dat(:);
    end
    fwrite(fidO,dat,'int16');
end

remainder = nBytes/(sizeInBytes*nbChan) - nbChunks*chunk;
if ~isempty(remainder)
    dat = fread(fidI,nbChan*remainder,'int16');
    if ~isempty(channelIx)
        dat = reshape(dat,[nbChan remainder]);
        dat = dat(channelIx,:);
        dat = dat(:);
    end
    fwrite(fidO,dat,'int16');
end
close(h);

fclose(fidI);
fclose(fidO);




