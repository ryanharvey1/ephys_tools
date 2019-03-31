function spike_count = SMatrix_getSpikeCounts(S)
% returns a vector of spike_counts (number of spikes in corresponding ts-object) of all cells in Spike-Matrix S
%
%  spike_count = SMatrix_getSpikeCounts(S)
%
%  S ... a Spike matrix (S-matrix) =  a cell array of ts-objects of spiketimes.
%  spike_count ... vector with same size as S holding the number of spikes of each ts-object in S.
%
% PL Feb. 2003

spike_count = zeros(size(S));
for i = 1:length(S)
    spike_count(i) = length(data(S{i}));        
end