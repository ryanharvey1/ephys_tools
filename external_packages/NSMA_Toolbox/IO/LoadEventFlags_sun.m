function [evtsd, ts, param1, param2, bytes] = LoadEventFlags_sun(fname)
%
% load a Cheetah_sun event flag file 
%  
% INPUT: fname ..... event flag filename string 
% 
% OUTPUT:  evtsd ... timestamped tsd object with Data(evtsd) =  TTL output values array
%
% PL Nov. 1999
% Version 0.1

RECSIZE = 20;     % Record Size (bytes) of a Cheetah NT event record

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
[tmpdn, tmpfn] = fileparts(fname);




%--- get EV file size
H = ReadHeader(fp);
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
ts = zeros(nRec, 1); 
bytes = zeros(nRec,20);

param1 = zeros(nRec,5);
param2 = zeros(nRec,3);
for i = 1:nRec 
   bytes(i,:) = fread(fp, [1, 20], 'char');
   % swap bytes 
    param1(i,1)  = bytes(i,1)*256 + bytes(i,2);
    param1(i,2)  = bytes(i,3)*256 + bytes(i,4);
    param1(i,3)  = bytes(i,5)*256 + bytes(i,6);
    param1(i,4)  = bytes(i,7)*256 + bytes(i,8);
    param1(i,5)  = bytes(i,9)*256 + bytes(i,10);
    ts(i)        = (bytes(i,11)*256 + bytes(i,12))*2^16 +  bytes(i,13)*256 + bytes(i,14);
    param2(i,1)  = bytes(i,15)*256 + bytes(i,16);
    param2(i,2)  = bytes(i,17)*256 + bytes(i,18);
    param2(i,3)  = bytes(i,19)*256 + bytes(i,20);
end%for

% unset bit 16 in param2(:,2)
param2(:,2) = bitset(param2(:,2),16,0);

evtsd = tsd(ts,param2(:,2));