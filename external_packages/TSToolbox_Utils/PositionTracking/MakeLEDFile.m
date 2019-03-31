function [returnVar,msg] = MakeLEDFile(fbasename,varargin)

% USAGE:
%     [returnVar,msg] = MakeLEDFile(filebasename,channel)
% INPUTS:
%     filebasename: file base name
%     channel (optionnal): sync pulse channel. Default is last channel in the eeg file.
% 
% Adrien Peyrache 2011

msg = '';
returnVar = 0;

try 

    info = xml_load([fbasename '.xml']);
    nChan = str2double(info.acquisitionSystem.nChannels);
    
    if ~isempty(varargin)    
        channel = varargin{1};
        if ~isnumeric(channel)
            error('Channel Nb should be numeric');
        end
    else
        channel = nChan;
    end

    fs = str2double(info.fieldPotentials.lfpSamplingRate);    
    led = LoadBinary([fbasename '.eeg'],'channels',channel,'nChannels',nChan,'frequency',fs);
    
    fid = fopen([fbasename '.led'],'w');
    fwrite(fid,led,'int16');
    fclose(fid);
    returnVar = 1;
    
catch
  msg = [msg '\n' lasterr];      
  returnVar = 0;
end
    
        
        