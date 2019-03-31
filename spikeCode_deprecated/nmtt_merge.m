% Script for parsing Neuralynx video (tracked) and tt files.
% dahamilt@unm.edu
% Created : 5/18/2015
% Last Modified : 5/18/2015
%
% DESCRIPTION : This script processes video (tracked) animal position/head
% direction data and single unit data (tt files) to create firing rate
% maps. A firing rate map is generated for each tt file. 
%
% USE : Navigate to the directory containing the data. 
% The script requires two arguments. The nvt filename, the spikesorted tetrode filename. The
% latter must have a .txt extension. The script currently processes the nvt
% file in its native binary format.

% e.g., 
% nmtt_merge('VT1.nvt', 'TT3_SS_01.txt')


function [ outputs ] = nmtt_merge(filename_vt, filename_ss)

vtd = ReadNVT(filename_vt); % video tracked data
ssfid = fopen(filename_ss, 'rt');
ssd = textscan(ssfid, '%u64'); %spike sorted data
fclose(ssfid);
ssd = ssd{1};

% assign ts, x, and y data from the video record
ts = vtd.TimeStamp;
ts = double(ts);
ssd = double(ssd);
xs = vtd.Xloc;
ys = vtd.Yloc;

[l,w] = size(ssd); % l = length of ssd array

ts(1)
ssd(1)

ts(1) - ssd(1)
for i = 1:l
    diffary = abs(ts - ssd(i)); % subtract the ss timestamp from all vt timestamps
    [m,ind] = min(diffary);    %find the min (closest timestamp) and its element
    diffs(i) = m;
    ssx(i) = xs(ind); % store x for spike
    ssy(i) = ys(ind); % store y for spike
end

ssx;

% get range for plot
minx = min(xs)
maxx = max(xs);
miny = max(ys);
maxy = max(ys);


figure(1);
scatter(ssx,ssy);
%axis([minx maxx miny maxy]);


min(diffs);
max(diffs);
mean(diffs);

