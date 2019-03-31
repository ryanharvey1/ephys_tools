function [T] = GenerateSyncedNeuralTimes(clkFname, dataFname)
% Generates a matlab-time synced time vector for the Intan neural data
% 
% USAGE
%
%   GenerateSyncedNeuralTimes(clkFname, dataFname)
%
% INPUT
%   clkFname - Intan DIN (digital input) file w/ clock signal
%   dataFname - Intan DIN (digital input) file w/ data signal
%
% OUTPUT
%   T - aligned time vector
%
% May 2017
% Brendan Hasz
% haszx010@umn.edu
% David Redish Lab, University of Minnesota Twin Cities

function N = bi2num(D)
    % Convert a binary vector to number
    N = sum(D.*power(2,0:length(D)-1));
end
    
% Get the CLOCK data file
fid = fopen(clkFname, 'rb');
sz = dir(clkFname);
CLK = fread(fid, sz.bytes/2, 'int16');
fclose(fid); %close the file

% Get the DATA data file
fid = fopen(dataFname, 'rb');
sz = dir(dataFname);
DATA = fread(fid, sz.bytes/2, 'int16');
fclose(fid);

%Get rising/falling clock edge inds
cinds = 1+find(diff(CLK)); 

% Get inds of each timestamp, and their binary values
Nts = length(cinds)/34; %number of timestamps
tsinds = zeros(Nts, 1); %index at the start of each timestamp
tsvals = zeros(Nts, 1); %values of each timestamp (ms since MatlabTime=0)
for ts = 1:Nts
    st_ind = (ts-1)*34+1; %index of the 1st pulse of this timestamp
    tsinds(ts) = cinds(st_ind); %index @ start of first pulse in this t.s.
    bvals = DATA(cinds(st_ind+2:st_ind+33)); %data value @ each clk edge for this t.s.
    tsvals(ts) = bi2num(bvals')/1000; %convert from binary to decimal, and convert from ms to seconds
end

% Generate time vector with time values corresponding to Matlab time and
%   indexes corresponding to the neural data
T = zeros(size(DATA));
T(1:tsinds(1)) = linspace(tsvals(1)-tsinds(1)/30000, tsvals(1), tsinds(1));
for ts = 2:Nts
    T((tsinds(ts-1)+1):tsinds(ts)) = linspace(tsvals(ts-1), tsvals(ts), tsinds(ts)-tsinds(ts-1));
end
T((tsinds(end)+1):end) = linspace(tsvals(end)+1/30000, tsvals(end)+(length(DATA)-tsinds(end))/30000, length(DATA)-tsinds(end));


end