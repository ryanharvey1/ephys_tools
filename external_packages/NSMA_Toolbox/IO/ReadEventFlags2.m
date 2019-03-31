function evt = ReadEventFlags2(fname)
% ReadEventFlags  Reads a Cheetah_NT event flag file and ouputs a struct (for people
% who don't like tsd's) with the more useful parts of the info in Cheetah event files (which contain a lot of air).
%
% evt = ReadEventFlags2(fname)
%
% INPUT: 
%   fname = event flag filename string 
% OUTPUTS:  
%   
%   evt = stuct of parallel column vectors (all same length) with fields:
%       evt.ts = event timestamps in NSMA timstamp units (0.1 msec units)
%       evt.port =  Cheetah port number on which event was received (on Cheetah 160 systems it is either 1 of 2 TTL ports (119 or ???) or the srting input port = 4) 
%       evt.code_dec = event codes in decimal numeric values
%       evt.code_hex = event codes in HEX string values
%       evt.code_bin =  event codes in binary string values
%       evt.code_dec_noStrobeBit = event codes in decimal value with the strobe bit
%           15 removed (since bit 15 is raised in ALL cheetah TTL event codes to trigger
%           the events)
%       evt.code_hex_noStrobeBit = event codes in decimal value with the strobe bit
%           removed (as above)
%       evt.string = cell array of strings of event strings (as given on the Cheetah
%           UserEvent Input Window or a string generated for each event received by a
%           TTL port).
%
% PL 05/2008, version 0.2


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

% repack relevant info into output struct
evt.ts = events_ts;
evt.port =  param2(:,1);        %Cheetah port number on which event was received (on Cheetah 160 systems it is either 1 of 2 TTL ports (119 or ???) or the srting input port = 4) 
evt.code_dec = param2(:,2);    % event codes in decimal numeric values
evt.code_hex = dec2hex(param2(:,2));       % event codes in HEX string values
evt.code_bin = dec2bin(param2(:,2));        %event codes in binary string values

% unset bit 16 in param2(:,2)
noStrobeCode = bitset(param2(:,2),16,0);

evt.code_dec_noStrobeBit = noStrobeCode; %event codes in decimal value with the strobe bit
%           15 removed (since bit 15 is raised in ALL cheetah TTL event codes to trigger
%           the events)
evt.code_hex_noStrobeBit = dec2hex(noStrobeCode); %event codes in decimal value with the strobe bit
%           removed (as above)
evt.code_bin_noStrobeBit = dec2bin(noStrobeCode); %event codes in binary string values with the strobe bit
%           removed (as above)
evt.string = evstring; %cell array of strings of event strings (as given on the Cheetah
%           UserEvent Input Window or a string generated for each event received by a
%           TTL port).

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
