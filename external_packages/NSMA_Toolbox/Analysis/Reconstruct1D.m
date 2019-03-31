function [rMu, rCI] = Reconstruct1D(S, SF, rStart, rEnd, rDT)

% [rMu, rCI] = Reconstruct1D(S, SF, rStart, rEnd, rDT)
%
% Does 1D reconstruction.
%
% INPUTS:
%     S = cell array of ts objects containing spike times
%     SF = cell array of scatter fields 
%     rStart, rEnd = start and end times within which to reconstruct
%     rDT = timeslices to take (measured in timestamps, i.e. 1000 ts = 100 ms = 0.1 sec)
%
% OUPUTS:
%     rMu = reconstructed position (tsd)
%     rCI = confidence interval on rMu (tsd)
%
% ALGORITHM

% ADR 1998
% version L4.0
% status UNDER CONSTRUCTION

StatusWarning('UNDER CONSTRUCTION', 'Reconstruct1D');

%------------------

% temporary, can we do with not binning position?
allBins = 1:25:250;
nBins = length(allBins);

%------------------

if (length(SF) ~= length(S))
   error('Size mismatch between SF & S.');
end

% take out empty cells
nCells = length(S);
keep = [];
for iC = 1:nCells
   if ~isempty(Data(S{iC})) & ~isempty(Data(SF{iC}))
      keep = cat(1, keep, iC);
   end
end
S = S(keep); SF = SF(keep);
nCells = length(S);

%-------------------     
% Calculate spike{c}(t) = 1 if spike occurs at time t, else 0

for iC = 1:nCells
   Spike{iC} = MakeQfromS(S(iC), rDT, 'T_start', rStart, 'T_end', rEnd, 'ProgressBar', 'none');
end

%-------------------
% Calculate ISI{c}(t) = time since last spike from time t
for iC = 1:nCells   
   SpikeTimes = Data(Restrict(S{iC}, rStart, rEnd));
   IDofLastSpike = cumsum(Data(Spike{iC}));
   if isempty(SpikeTimes) | Data(Spike{iC},rStart+1) == 0
      SpikeTimes = cat(1, Data(S{iC}, rStart), SpikeTimes);
      IDofLastSpike = 1 + IDofLastSpike;
   end   
   TimeofLastSpike{iC} = SpikeTimes(IDofLastSpike);
end
return;
%=============================================================
muC = zeros(nCells,1);            % mean for each cell's field
sigmaC = zeros(nCells,1);         % sigma for each cell's field
varC = zeros(nCells, 1);          % variance for each cell's field (speeds computation to precomput)

% find gaussian parameters from scatter fields
for iC = 1:nCells
   [muC(iC), sigmaC(iC)] = normfit(Data(SF{iC}));   
end
varC = sigmaC .* sigmaC;

% find poissfit parameters from scatter fields
for iC = 1:nCells
   lambdaC(iC) = poissfit(diff(Range(SF{iC}, 'ts')));
end

% calculate Q matrix
Q = MakeQfromS(S, rDT, 'T_start', rStart, 'T_end', rEnd, 'ProgressBar', 'none');
allTime = Range(Q,'ts');
nTime = length(allTime);

% assume uniform occupancy
logPX = 1/nBins;

% generate data for each timestep
logPS = zeros(nCells,1); logPSX = zeros(nCells,1); 
for iT = 1:nTime
   DisplayProgress(iT,nTime);
   DQ = Data(Q,allTime(iT));
   for iX = 1:nBins
      for iC = 1:nCells
         DisplayProgress(iC,nCells,'UseGraphics',0);
         lambda = normpdf(allBins(iX)-muC(iC),muC(iC),sigmaC(iC));
         if DQ(iC)
            logPSX(iC) = logpoisspdf(allTime(iT)-Data(S{iC},allTime(iT)),1/(lambda+1e-100));
            logPS(iC) = logpoisspdf(allTime(iT)-Data(S{iC},allTime(iT)),lambdaC(iC));
         else
            logPSX(iC) = 1 - logpoisspdf(allTime(iT)-Data(S{iC},allTime(iT)),1/(lambda+1e-100));
            logPS(iC) = logpoisspdf(allTime(iT)-Data(S{iC},allTime(iT)),lambdaC(iC));   
         end % if spike
      end % for all cells
   end % for all bins
   logPXS = sum(logPSX) + logPX - logPS;
   
   % we need to fill this out ...
   PXS = SampleDistribution(exp(logPXS));
   
   % calculate rMu and rCI for each timestep
   [rMu(iT), rSigma, rCI(iT)] = normfit(PXS);
end % for all time
