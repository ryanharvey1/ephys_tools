% [Avs StdErr] = TriggeredAv(Trace, nBefore, nAfter, T, G)
%
% computes triggered averages from Trace at the times given by T.
% Trace may be 2D, in which case columns are averaged separately.
% The output will then be of the form Avs(Time, Column)
% nBefore and nAfter give the number of samples before and after
% to use.
%
% G is a group label for the trigger points T.  In this case the
% output will be Avs(Time, Column, Group)
%
% StdErr gives standard error in the same arrangement

function [Avs, StdErr] = TriggeredAv(Trace, nBefore, nAfter, T, G)

if (nargin<5 | length(G) == 1)
	G = ones(length(T), 1);
end

nColumns = size(Trace,2);
nSamples = nBefore + nAfter + 1;
nGroups = max(G);
maxTime = size(Trace, 1);

Avs = zeros(nSamples, nColumns, nGroups);

BlockSize = floor(200000/nSamples); % memory saving parameter

for grp = 1:nGroups

	Sum = zeros(nSamples, nColumns);
	SumSq = zeros(nSamples, nColumns);
	MyTriggers = find(G==grp & T > nBefore & T <= maxTime-nAfter);
	nTriggers = length(MyTriggers);
	
	% go through triggers in groups of BlockSize to save memory
	for Block = 1:ceil(nTriggers/BlockSize)
		BlockTriggers = MyTriggers(1+(Block-1)*BlockSize:min(Block*BlockSize,nTriggers));
		nBlockTriggers = length(BlockTriggers);
		
		TimeOffsets = repmat(-nBefore:nAfter, nBlockTriggers, 1);
		TimeCenters = repmat(T(BlockTriggers), 1, nSamples);
		TimeIndex = TimeOffsets + TimeCenters;
		
		Waves = Trace(TimeIndex,:);
		Waves = reshape(Waves, [nBlockTriggers, nSamples, nColumns]);
		Sum = Sum + reshape(sum(Waves, 1), [nSamples, nColumns]);
		SumSq = SumSq + reshape(sum(Waves.^2,1), [nSamples, nColumns]);
	end
	
	Avs(:,:,grp) = Sum/nTriggers;
	StdErr(:,:,grp) = sqrt(SumSq/nTriggers - Avs(:,:,grp).^2) / sqrt(nTriggers);
end