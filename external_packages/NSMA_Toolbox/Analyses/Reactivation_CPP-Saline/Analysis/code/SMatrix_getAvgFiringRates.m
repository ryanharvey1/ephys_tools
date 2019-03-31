function avg_rate = SMatrix_getAvgFiringRates(S)
% returns a vector of average firing rates (number of spikes in corresponding ts-object divided by length of S epoch) of all cells in Spike-Matrix S
%
%  avg_rate = SMatrix_getAvgFiringRates(S)
%
%  S ... a Spike matrix (S-matrix) =  a cell array of ts-objects of spiketimes.
%  avg_rate ... vector with same size as S holding the number of spikes of each ts-object in S divided 
%               by the length (in sec) of the S-epoch as returned by [SessStartTS, SessEndTS] = SessionStartEndTS(S).
%               There is no check for gaps in the recording session. The units of avg_rate are in Hz. 
%
% PL Feb. 2003

% get epoch length (sec)
[SessStartTS, SessEndTS] = SessionStartEndTS(S);
interval_sec = (SessEndTS - SessStartTS)/10000;

% get spike counts and compute rate
avg_rate = zeros(size(S));
for i = 1:length(S)
    avg_rate(i) = length(data(S{i}))/interval_sec;        
end
