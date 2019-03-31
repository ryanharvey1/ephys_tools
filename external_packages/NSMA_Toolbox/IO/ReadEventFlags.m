function [evtsd, events_ts, param1, param2, param3, evstring] = ReadEventFlags(fname)

% ReadEventFlags  Reads a Cheetah_NT event flag file
%
% [evtsd, events_ts, param1, param2, param3, evstring] = ReadEventFlags(fname)
%
% INPUT: 
%   fname = event flag filename string 
% OUTPUTS:  
%   evtsd = timestamped tsd object with Data(evtsd)=TTL output values array
%   events_ts = event timestamps (same as in evtsd)
%   param1,2,3 = components of events file - evtsd data comes from param2
%   evstring = cell array of strings of events names
%
% PL 11/99, version 0.1, last modified 11/99 by PL


RECSIZE = 184;     % Record Size (bytes) of a Cheetah NT event record

%---- Open EV file
if ~exist(fname, 'file')
   errordlg(['No Event Flag file "' fname '" found! Check path and try again!!!'], ...
         'LoadEventFlags Error');
   return;    %abort
end%if
[fp, msg] = fopen(fname,'rb','n');
if (fp == -1) 
   error(['Error opening file "' fname '" : ' msg]); 
end%if
header = SkipCheetahNTHeader(fp);
if(length(header)>= 2000)
    fprintf('%s',header(1:2000));
end
[tmpdn, tmpfn] = fileparts(fname);


%--- get EV file size
fStart = ftell(fp);    % memorize start position after header
whambang = fseek(fp,0,'eof');   
if whambang 
   [mess, errnum] = ferror(fp); 
   errordlg(['Something is fatally wrong with current Eventflag file:' ...
         mess , ' ErrNo: ', num2str(errnum)], 'LoadEventFlags Error'); 
   return; 
end%if 
eofPos = ftell(fp);
fSize = eofPos - fStart;    % file size in bytes
fseek(fp,fStart,'bof');     % rewind to file start position (after header)
nRec = floor(fSize/RECSIZE);
lastRecSize = mod(fSize,RECSIZE);
if lastRecSize > 0
   warning(sprintf('File %s has an incomplete last record of size %i bytes. It is ignored!',...
      tmpfn, lastRecSize));
end%if
disp(sprintf('Loading file %s with %i event records of size %i bytes.',...
      tmpfn, nRec, RECSIZE));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  main loop over records
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
events_ts = zeros(nRec, 1); 
param1 = zeros(nRec,3);
param2 = zeros(nRec,5);
param3 = zeros(nRec,8);
evstring = cell(nRec,1);    
for i = 1:nRec 
   param1(i,:)  = fread(fp, [1,3] , 'int16' ); 
   events_ts(i) = fread(fp,     1 , 'int64' )/100;   
   param2(i,:)  = fread(fp, [1,5] , 'uint16' );
   param3(i,:)  = fread(fp, [1,8] , 'uint32' );
   evstring{i}  = char(fread(fp, [1,128] , 'char' ));
end%for

% unset bit 16 in param2(:,2)
param2(:,2) = bitset(param2(:,2),16,0);

evtsd = tsd(events_ts,param2);

fclose(fp);

return



%% -----------------------------------------------------------

function headercontent = SkipCheetahNTHeader(fp)
%
% Cheetah NT versions 1.3 or greater have a standard Neuralynx header of ascii strings of total
% length of 16384 bytes in the beginning of the binary file. 
% Check if Cheetah header present (in versions > 1.3) and skip header
%
% returns headercontent string if new Neuralynx header is present and was skipped
%    or   empty string '' if NO  header is present

headercontent = '';      %% if first 8 bytes of AnyCheetahfile.dat contain "########" there is a header 
                         %% of length 16,384 bytes to skip (from the very beginning including the eight '#') 
NHEADERBYTES = 16384;	

headerflag = char(fread(fp, [1,8], 'uchar'));  
fprintf( 'Headerflag =  %s \n', headerflag);
fseek(fp,0,'bof');       %% reset filepointer to the beginning of file
if (strcmpi(headerflag,'########') == 1)
   headercontent = char(fread(fp,[1,NHEADERBYTES],'uchar')); 
   %%fseek(fp,NHEADERBYTES,'bof');  %% set filepointer after byte NHEADERBYTES
   fprintf('NT-Header read (%d bytes)\n',NHEADERBYTES);
   return;
end 
return
