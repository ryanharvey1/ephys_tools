function WV = ExtractWaveforms(TT, S, varargin)

% WV = ExtractWaveforms(TT, S)
%
% INPUTS:
%     TT = tsd struct(t = timestamps of all spikes, data = waveforms of all spikes)
%     S = ts structure (or cell array of same)
% 
% OUTPUTS:
%     WV = tsd containing only those spikes in each cell S

% ADR 1998
% version L4.0
% status PROMOTED

nTrodes = 4;
nSamples = 32;
Extract_varargin;

switch class(S)
case 'cell'
   WV = {};
   for iC = 1:length(S)
      WV{iC} = Restrict(TT, Data(S{iC}));
   end
case 'ts'   
   WV = Restrict(TT, Data(S));
otherwise
   error('Unknown input type.');
end
