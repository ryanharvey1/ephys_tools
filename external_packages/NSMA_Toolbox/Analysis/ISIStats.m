function [mISI, vISI, nSp] = ISIStats(S)

% [mISI, vISI, nSp] = ISIStats(S)
% 
% INPUTS
%     S = cell array of ts objects
%
% OUTPUTS
%     mISI = mean interspike interval
%     vISI = variance of interspike interval
%     nSp = number of spikes fired

% ADR 1998
% version vL4.0
% status PROMOTED

if ~isa(S, 'cell')
   error('Type Error: S is not of type cell.');
end

nCells = length(S);
mISI = zeros(size(S));
vISI = zeros(size(S));
nSp = zeros(size(S));

for iC = 1:nCells
   if ~isa(S{iC}, 'ts')
      error(['Type Error: S{' num2str(iC) '} is not a ts object.']);
   end   
   mISI(iC) = mean(diff(Data(S{iC})));
   vISI(iC) = var(diff(Data(S{iC})));
   nSp(iC) = length(Data(S{iC}));
end
