function sparsity = Sparsity(varargin)
% Computes the sparsity of a cell
%
%     Calculate single-cell sparsity index (equivalently
%        'lifetime sparseness') as defined in Ahmed and Mehta (Trends in
%        Neuroscience, 2009)
%
%   S=(1/n sum(ri))^2 / (1/n sum(ri^2))
%
% Where i is a spatial bin and ri is the firing rate of the cell in bin i of an environment
% containing a total or n spatial bins. A sparsity value of 1 implies no single-cell
% sparseness. A sparsity value approaching 0 is indicative of maximal single-cell
% sparseness, and implies a greater amount of spatial information in each spike emitted by
% that cell.
%
% See parameters below.
%
%
%   PARAMETERS
%
%   occupation_thresh   Bins that had less than this number of seconds of
%                       occupation are not included in the score. (0)
%   ratemap             pre-computed ratemap
%
%   occupancy           occupancy map
%
% Ryan Harvey 2019

p = inputParser;

p.addParameter('occupation_thresh', 0, @isnumeric);
p.addParameter('ratemap', [], @isnumeric)
p.addParameter('occupancy', [], @isnumeric)

p.parse(varargin{:});

occupation_thresh = p.Results.occupation_thresh;
ratemap = p.Results.ratemap;
occupancy = p.Results.occupancy;

% normalize to probability
occupancy = occupancy / sum(occupancy(:));

ratemap(occupancy <= occupation_thresh) = NaN;

ratemap=ratemap(:);

num = (nansum(ratemap)./ length(ratemap)).^2;
den = nansum(ratemap.^2)./ length(ratemap);
sparsity = num / den;
end




