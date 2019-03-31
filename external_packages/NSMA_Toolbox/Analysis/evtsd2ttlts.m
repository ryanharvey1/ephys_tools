function [ttl, ttlhex, A_hex, B_hex, A_dec, B_dec] = evtsd2ttlts(evtsd) 
%
% [ttl, ttlhex, A_hex, B_hex, A_dec, B_dec] = evtsd2ttlts(evtsd) 
% [ttlts, ttlhex] = evtsd2ttlts(evtsd)
%
% convert eventflags tsd into human readable cell array ttlts
%
% INPUT: 
%     evtsd ... eventflag tsd (from LoadEventFlags) 
% OUTPUT:
%     ttls ...  {nevents x 2} cellarray: { timstamp  hex-ttl-codes }
%     ttlhex ... cellarray of hex-ttl-codes strings
%
%  PL Nov 1999
%  version 0.1

%unpack evtsd
ts = Range(evtsd,'ts');
param2 = Data(evtsd);
ttldec = param2(:,2);

% split ttl word into low (A) and high (B) bytes
byteAmask = hex2dec('00FF');
byteBmask = hex2dec('FF00');
A_dec = bitand(ttldec,byteAmask);
B_dec = bitshift(bitand(ttldec,byteBmask),-8);

A_hex = cellstr(dec2hex(A_dec));
B_hex = cellstr(dec2hex(B_dec));
ttlhex = cellstr(dec2hex(ttldec));
ttl = [num2cell(ts) ttlhex B_hex num2cell(A_dec)];