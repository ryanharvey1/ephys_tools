function [data, norm, skip, raw]  = ReadNicoletScopeFile(fname)

% ReadScopeFile  Reads propriety binary file from Nicolet 310 Storage Oscilloscope
%
% [data, norm, skip, raw]  = ReadScopeFile(fname)
%
% INPUTS:
%       fname - matlab string with the filename
% OUTPUTS:
%       data - 4000x2 data matrix:  column 1 = time (msec), column 2 = Amplitude (mVolts)
%       norm - ? <-(another programmer who comes along
%       skip - ? <- later who knows what these are can 
%       raw - ?  <- help fill in these descriptions)
%
% Reads a propriety binary file from the Nicolet 310 Storage Oscilloscope
% and converts and returns the data in a 2 column array 'data'


fclose('all');
fname
[fid,message] = fopen(fname,'r');
disp(message);

[norm, count] = fread(fid, [5,20], 'uchar');
%fprintf('count = %d\n', count);
norm = char(norm');
norm = str2num(norm);
vzero = norm(3);
hzero = norm(4);
vnorm = norm(5)*10^norm(6);
hnorm = norm(7)*10^norm(8)*1000;

[skip, count] = fread(fid, [1,92], 'uchar');
%fprintf('count = %d\n', count);

[raw, count] = fread(fid, 4000, 'int16');
%fprintf('count = %d\n', count);

data(:,1) = ((0:length(raw)-1)'-hzero) * hnorm;
data(:,2) = (raw-vzero) * vnorm;

fclose(fid);




