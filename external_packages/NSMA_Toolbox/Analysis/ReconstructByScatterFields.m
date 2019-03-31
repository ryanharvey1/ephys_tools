function [rMu, rSigma, xbins, R] = ReconstructByScatterFields(S, SF, WT, V, rStart, rEnd, rDT, varargin)

% [rMu, rSigma, xbins, R] = ReconstructByScatterFields(S, SF, WT, V, rStart, rEnd, rDT, parameters)
%
% INPUTS
%    S = cell array of ts objects (spike times)
%    SF = cell array of scatter fields
%    WT = tsd of dimension along which scatter fields exist (dimension with which to reconstruct)
%    rStart = timestamp at which to start reconstruction
%    rEnd = timestamp at which to end reconstruction
%    rDT = timestep for reconstruction
%
% OUTPUTS
%    rMu = ctsd of mean reconstructed position at each timestep
%    rSigma = ctsd of mean reconstructed position at each timestep
%    xbins = bins for R
%    R = ctsd of position reconstruction binned output
%    note: if function is called with 1 or 2 outputs, then R is not generated
%          which will speed the function up considerably
%
% PARAMETERS
%       sigma, stddev for guassian weighting function around t (default 1000 ts = 100 ms)
%       NormalizeByTotalSpikes (default 1) 
%       NormalizeByOccupancy (default 1)
%       xBins bins to split R into (default linspace(min, max, 100))
%       Radius, radius to normalize occ by (default 10)
%
% ALGO
%    For each time slice t, fits a gaussian to the distribution 
%        integral_0^t (sum_c (SF_c(tau) * WT_c * S_c(tau)) * alpha(t - tau) dt)
%    where t = time of time slice
%          c = cell index
%          SF = scatter field
%          WT = weight of scatter field
%          S_c = 1 if spike at time tau, 0 otherwise
%          alpha = gaussian smoothing function with stddev sigma
%
% ADR 1998
% version L5.0 
% status STILL TESTING

% L4.1 now can use weighted scatter fields
% L4.1 now uses a gaussian around t instead of exponential leading up to t
% L5.0 no longer uses weighted scatter fields, does the normalization here.

StatusWarning('STILL TESTING', 'ReconstructByScatterFields');

%------------------------
% parameters

adrlib;
sigma = 1000; 
xBins = linspace(min(Data(V)), max(Data(V)), 100);
NormalizeByTotalSpikes = 1;
NormalizeByOccupancy = 1;
Radius = 10;
Extract_varargin;

if nargout == 4
   flagReturnFullReconstruction = 1;
   nBins = length(xBins);
else
   flagReturnFullReconstruction = 0;
end

nCells = length(S);
timeSlices = rStart:rDT:rEnd;
nTimeSlices = length(timeSlices);

alphaToZero = 3 * sigma;  % alpha is close to zero here, we can ignore it

%------------------------
% first build a huge 'tsd' of when each cell fires

spikeTimes = []; cellIDs = []; 
for iC = 1:nCells
   DisplayProgress(iC, nCells, 'Title', 'Preprocessing');
   S0 = Restrict(S{iC}, rStart - alphaToZero, rEnd + alphaToZero);
   spikeTimes = cat(1, spikeTimes, Data(S0));
   cellIDs    = cat(1, cellIDs,    repmat(iC, size(Data(S0))));
end
[spikeTimes, order] = sort(spikeTimes);
cellIDs = cellIDs(order);
clear order S0

%-------------------------
% next for each time slice, calculate eqn
% carry equation out to alphaToZero

if flagReturnFullReconstruction
   R = zeros(nTimeSlices, nBins);
end

for iT = 1:nTimeSlices
   
   DisplayProgress(iT, nTimeSlices);
   
   f = find((spikeTimes > timeSlices(iT)-alphaToZero) & (spikeTimes < timeSlices(iT)));
   slicedSpikeTimes = spikeTimes(f);
   slicedCellIDs = cellIDs(f);
   x = [];, w = [];
   for iS = 1:length(slicedSpikeTimes)
      x = cat(1, x, DATA(SF{slicedCellIDs(iS)}));
      alpha = normpdf(slicedSpikeTimes(iS), timeSlices(iT), sigma);
      w = cat(1, w, alpha * DATA(WT{slicedCellIDs(iS)}));
   end
   
   if isempty(x)
      rMu(iT) = NaN;
      rSigma(iT) = NaN;
   else      
      rMu(iT) = nansum(w .* x)/nansum(w);
      rSigma(iT) = sqrt(nansum(w .* (x - rMu(iT)).^2)/nansum(w));
      
      if flagReturnFullReconstruction
         for ik = 2:nBins
            f = find(x > xBins(ik-1) & x <= xBins(ik));
            if ~isempty(f)
               R(iT, ik) = w(f)' * x(f);
            end
         end
      end
      
   end
   
end

%-----------------------
% fill ctsd's

rMu = ctsd(rStart, rDT, rMu');
rSigma = ctsd(rStart, rDT, rSigma');

if flagReturnFullReconstruction
   R = ctsd(rStart, rDT, R);
end