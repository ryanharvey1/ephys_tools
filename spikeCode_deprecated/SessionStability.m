function [rSpatial, vFirst, vSec] = SessionStability(data,bins,sampleRate)
% stability
%
% Input: 
%       - smoothed matrix
%
% Output:
%       - r = pearsons correlation coefficient
% 
% 
% Created by Ben C June 2014; requires Chronux


% extract x y coords and smooth
rxf = data(:,1);
ryf = data(:,2);
rxs = runline(rxf,3,1); % smooth with 3pt window and 1pt step (based on Valente et al 2007 PloS One)
rys = runline(ryf,3,1);

% extract spike locations
spks = data(:,5);

% determine start and end of first and second half of session
timestamps = 1:length(data(:,1));
sessLength_samples = numel(timestamps);
First_End = round(sessLength_samples/2);
Sec_Start = First_End + 1;
Sec_End = sessLength_samples;

% determine spikes and path in first and second half
First_spks = spks(1:First_End);
Sec_spks = spks(Sec_Start:Sec_End);
First_rxs = rxs(1:First_End);
First_rys = rys(1:First_End);
Sec_rxs = rxs(Sec_Start:Sec_End);
Sec_rys = rys(Sec_Start:Sec_End);
First_spksi = find(First_spks == 1);
First_spks_xs = rxs(First_spksi);
First_spks_ys = rys(First_spksi);
Sec_spksi = find(Sec_spks == 1);
Sec_spks_xs = rxs(Sec_spksi);
Sec_spks_ys = rys(Sec_spksi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bin the positional data for the first half
MinY = min(First_rys);
MaxY = max(First_rys);

MinX = min(First_rxs);
MaxX = max(First_rxs);

edges{1} = linspace(MinY, MaxY, bins+1);
edges{2} = linspace(MinX, MaxX, bins+1);

% bin positional data
occMatrix1 = [First_rys, First_rxs];
Omatrix1 = hist3(occMatrix1,'Edges',edges);
Omatrix1(1,:) = [];
Omatrix1(:,end) = [];
occ1 = Omatrix1/sampleRate;

% bin spike data
spikeMatrix1 = [First_spks_ys, First_spks_xs];
Smatrix1 = hist3(spikeMatrix1,'Edges',edges);
Smatrix1(1,:) = [];
Smatrix1(:,end) = [];

% divide binned spikes by occupancy to get rate maps and store data
BFiringRateMatrix1 = Smatrix1./occ1;
FilledRateMatrix1 = BFiringRateMatrix1;
FilledRateMatrix1(isnan(FilledRateMatrix1)) = 0;
FilledRateMatrix1(isinf(FilledRateMatrix1)) = 0;

% smooth raw rate histogram for first half of session
npts = 1;
gsize = 5;
[vFirst] = smooth(FilledRateMatrix1, npts, gsize, npts, gsize);
% vFirst_NaN = 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bin the positional data for the second half
occMatrix2 = [Sec_rys, Sec_rxs];
Omatrix2 = hist3(occMatrix2,'Edges',edges);
Omatrix2(1,:) = [];
Omatrix2(:,end) = [];
occ2 = Omatrix2/sampleRate;

% bin spike data
spikeMatrix2 = [Sec_spks_ys, Sec_spks_xs];
Smatrix2 = hist3(spikeMatrix2,'Edges',edges);
Smatrix2(1,:) = [];
Smatrix2(:,end) = [];

% divide binned spikes by occupancy to get rate maps and store data
BFiringRateMatrix2 = Smatrix2./occ2;
FilledRateMatrix2 = BFiringRateMatrix2;
FilledRateMatrix2(isnan(FilledRateMatrix2)) = 0;
FilledRateMatrix2(isinf(FilledRateMatrix2)) = 0;

% smooth raw rate histogram for second half of session
npts = 1;
gsize = 5;
[vSec] = smooth(FilledRateMatrix2, npts, gsize, npts, gsize);

% reshape data so it is column vector
% ai = reshape(FilledRateMatrix1,bins*bins,1);
% bi = reshape(FilledRateMatrix2,bins*bins,1);

% % remove bins with zeros
% zero_ai = find(ai < 0);
% ai(zero_ai,:) = NaN;
% bi(zero_ai,:) = NaN;
% zero_bi = find(bi < 0);
% bi(zero_bi,:) = NaN;
% ai(zero_bi,:) = NaN;

% obtain pearsons correlation between sessions
% ai = reshape(ai,bins,bins);
% bi = reshape(bi,bins,bins);
rSpatial = corr2(BFiringRateMatrix1,BFiringRateMatrix2);