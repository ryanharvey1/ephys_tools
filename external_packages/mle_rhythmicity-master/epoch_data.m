function [ data, inds ] = epoch_data( spk_ts, epochs, varargin )
%EPOCH_DATA Apply epochs to spike times for mle_rhythmicity
% 
% INPUT
%   spk_ts - Times of each spike
%   epochs - A mX2 matrix with the first column containing epoch starts and
%       the second containing epoch ends. 
%
% PARAMETERS
%   max_lag (0.6): Examination window
%   epochmode ('lead') - Can be 'lead' or 'lag'. If lead, excludes leading
%       spikes which would have lags falling out of the epoch. If lag, only
%       excludes lags falling out of epoch. 
%
% OUTPUT
%   data - Cell array of epoched lags for use in mle_rhythmcity
%   inds - Cell array of indicies of spikes at each lag for further
%       analysis
%
% See also mle_rhythmicity, rhythmicity_covar, rhythmicity_pdf
%
% Copyright 2015-2016 Trustees of Boston University
% All rights reserved.
%
% This file is part of mle_rhythmicity revision 2.0. The last committed
% version of the previous revision is the SHA starting with 93862ac...
%
% This code has been freely distributed by the authors under the BSD
% license (http://opensource.org/licenses/BSD2-Clause). If used or
% modified, we would appreciate if you cited our papers:
%
% Climer JR, DiTullio R, Newman EL, Hasselmo ME, Eden UT. (2014),
% Examination of rhythmicity of extracellularly recorded neurons in the
% entorhinal cortex. Hippocampus, 25:460-473. doi: 10.1002/hipo.22383.
%
% Hinman et al., Multiple Running Speed Signals in Medial Entorhinal
% Cortex, Neuron (2016). http://dx.doi.org/10.1016/j.neuron.2016.06.027

% Validate epochs
if size(epochs,2)~=2||any(epochs(:,1)>epochs(:,2))||any(arrayfun(@(t)any(t>epochs(:,1)&t<epochs(:,2)),epochs(:,1)))||any(arrayfun(@(t)any(t>epochs(:,1)&t<epochs(:,2)),epochs(:,2)))
    error('epoch_data:badEpochs','Error. Epochs must be a matrix of size mX2 as a list of non-overlapping epoch starts and stops.');
end

% Import data
ip = inputParser;
ip.addParamValue('epochmode','lead');
ip.addParamValue('max_lag', 0.6);
ip.parse(varargin{:});
for j = fields(ip.Results)'
    eval([j{1} ' = ip.Results.' j{1} ';']);
end

% Format lead_epochs based on epochmode
if isequal(lower(epochmode),'lead')
    lead_epochs = [epochs(:,1) epochs(:,2)-max_lag];
elseif iequal(lower(epochmode),'lag')
    lead_epochs = epochs;
else
    error('epoch_data:unrecognizedMode','Error. Unrecognized epochmode %s. Must be ''lead'' or ''lag''',epochmode);
end
lead_epochs = epochs(lead_epochs(:,2)>lead_epochs(:,1),:);% Keep only valid epochs

% Find the leading spikes
lead_spikes = spk_ts(arrayfun(@(t)any(t>=lead_epochs(:,1)&t<=lead_epochs(:,2)),spk_ts));

good_lags = arrayfun(@(t)any(t>=epochs(:,1)&t<=epochs(:,2)),spk_ts);

data = cell(size(lead_spikes));
inds = data;

for i=1:numel(data)
   data{i} = spk_ts(lead_spikes(i)<spk_ts&spk_ts-lead_spikes(i)<max_lag&good_lags)-lead_spikes(i); 
   inds{i} = find(lead_spikes(i)<spk_ts&spk_ts-lead_spikes(i)<max_lag&good_lags);
end

end

