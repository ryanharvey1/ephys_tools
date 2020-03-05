function mono_res = mono_synaptic(data) 
% mono_synaptic: get monosynaptic connections based on spike times
% ephys_tools wrapper for ce_MonoSynConvClick bz_PlotMonoSyn
% dependencies: Cell-Explorer (https://github.com/petersenpeter/Cell-Explorer)
%
% Ryan H 2020

[spikestimes,spikeIDs] = get_ce_MonoSynConvClick_input(data);


mono_res = ce_MonoSynConvClick(spikeIDs,spikestimes);

% mono_res = bz_PlotMonoSyn(mono_res);

end

function[spikestimes,spikeIDs] = get_ce_MonoSynConvClick_input(data)
spikestimes = vertcat(data.Spikes{:});
i = 0;
for s = 1:length(data.Spikes)
    idx = i+1:length(data.Spikes{s})+i;
    i = i+length(data.Spikes{s});
    spikeIDs(idx,1) =...
        repmat(str2double(extractBetween(data.spikesID.TetrodeNum(s),'TT','.mat')),...
        length(data.Spikes{s}),1);
    spikeIDs(idx,2) = repmat(data.spikesID.CellNum(s),length(data.Spikes{s}),1);
    spikeIDs(idx,3) = repmat(s,length(data.Spikes{s}),1);
end
end