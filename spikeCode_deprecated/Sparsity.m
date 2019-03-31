function [A,Ab,L] = Sparsity(R)

% R ... row cell array of rate histograms with nBins bins (column vectors) (weighted by occupancy)
% A ... Sparsity per X-bin over Population Vector, 
%        a la Treves and Rolls (numeric row vector indexed by icol=1:nBins)
%        (a measure of How man Cells are active at place i)
% Ab ... binary (left-right summed) version of Sparsity Array 
% L =   Sparseness of PlaceField, a measure for Localization (Concentration) of place field
%      L = numeric column vector indexed by cell number

% PL 1999

if ~isa(R,'cell')
   R={R};
end

% convert {1 x nCells} cell array of (nBins x 1) column vectors into (nBins x nCells) matrix
RC = cat(2,R{:});

[nBins, nCells] = size(RC);

A = sum(RC,2).^2./(nCells*sum(RC.^2,2));
L = sum(RC,1).^2./(nBins*sum(RC.^2,1));

mid = floor(nBins/2);
left_range = 1:mid;
right_range = mid+1:nBins;
RCb = [sum(RC(left_range,:),1) ; sum(RC(right_range,:),1)];
Ab = sum(RCb,2).^2./(nCells*sum(RCb.^2,2));

