function spatial_information = SpatialInformation(varargin)
% Computes the spatial information of a cell (bits/spike)
%
%
% See parameters below.
%
% See Cacucci et al 2007 Methods
%
%  PARAMETERS
%
%   occupation_thresh   Bins that had less than this number of seconds of
%                       occupation are not included in the score. (0)
%
%   n_thresh            min number of spikes before value is meaningless
%
%   ratemap             pre-computed ratemap
%
%   occupancy           occupancy map
%
%   n_spikes
% andrew 14 mat 2010
% enewman 20131004 added pre-computed ratemap parameter
% wchapman 20140514 recalculated 'F' by overall firing rate
% Modified by Ryan Harvey 2019

p = inputParser;

p.addParameter('occupation_thresh', 0, @isnumeric);
p.addParameter('n_thresh', 50, @isnumeric);
p.addParameter('ratemap', [], @isnumeric)
p.addParameter('occupancy', [], @isnumeric)
p.addParameter('n_spikes', [], @isnumeric)

p.parse(varargin{:});

occupation_thresh = p.Results.occupation_thresh;
ratemap = p.Results.ratemap;
n_thresh = p.Results.n_thresh;
occupancy = p.Results.occupancy;
n_spikes = p.Results.n_spikes;


occupancy = occupancy / sum(occupancy(:)); % normalize to probability

ratemap(occupancy <= occupation_thresh) = NaN;

F = nanmean(ratemap(:));

spatial_information = nansum(nansum(occupancy.*(ratemap./F).*log2(ratemap./F),1),2)./F;

spatial_information(n_spikes<n_thresh) = NaN; % where spiking was too low, get rid of spatial information score

end




