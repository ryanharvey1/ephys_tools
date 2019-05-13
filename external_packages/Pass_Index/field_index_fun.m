function [ field_index ] = field_index_fun( varargin )
% FIELD_INDEX_FUN - Calculated the field index.
%
% Calculates the field index along a trajectory. This function uses the
% pass_index_parser to generate its input structure.
%
% FIELD_INDEX = FIELD_INDEX_FUN(POS_TS,POS,SPK_TS)
% FIELD_INDEX = FIELD_INDEX_FUN(POS_TS,POS,SPK_TS,LFP_TS,LFP_SIG)
% FIELD_INDEX = FIELD_INDEX_FUN(POS_TS,POS,SPK_TS,LFP_TS,LFP_SIG,PARAMS)
% FIELD_INDEX = FIELD_INDEX_FUN(POS_TS,POS,SPK_TS,[],[],PARAMS)
%
%   ARGUMENTS
%   * POS_TS: Vector of time stamps for the sample state
%   * POS: MXN matrix of the sample state, where M is the number of samples
%   and N is the dimensions of POS
%   * SPK_TS: Spike times for the cell
%
%   OPTIONAL ARGUMENTS
%   * LFP_TS: Time stamps for the local field potential (LFP) Sample
%   * LFP_SIG: The LFP signal
%
%   PARAMETERS
%   * get_map (false): If true, returns a map instead of sampling
%   along the trajectory
%
%   PASS_INDEX_PARSER PARAMETERS
%   * method: Default 'grid'. Can be 'grid','place', or custom. Updates
%   other unset fields for these techniques.
%   * binside: Default 2*N, where N is the dimensionality of POS. Side of
%   the bins for rate mapping.
%   * smth_width: Default 3*BINSIDE, width of Gaussian smoothing kernel
%   * field_index: Default @field_index_fun, can be a vector of the same
%   number of elements as pos_ts, or can be a function handle which takes
%   in the same parameters as pass_index.
%   * sample_along: Default 'auto', can be 'arc_length', 'raw_ts', or a
%   nX2 matrix where n is the number of resampled steps, the first column
%   is the resampled timestamps and the second column is the sampled values,
%   or a function handle that returns a nX2 matrix as described above. Set
%   from 'auto' to 'arc_length' if method is 'place' or 'grid'.
%   * filter_band: Default 'auto', can be any positive frequency range in
%   cycles/unit sampled along using the ‘filter_band’ parameter.
%   Additionally, filter_band can be a function handle which returns a
%   modified signal. Set from 'auto' to [0.0749 0.0029] if 'method' is
%   'grid' and to [3*D 1/6*D].^-1, where D is the field width
%   determined by finding the N-dimensional volume of the region with at
%   least 10% of the maximum firing rate, and calculating the diameter of
%   the n-ball with the same volume.
%   * lfp_filter: Default [6 10]. can be changed to any frequency range in
%   Hz as [low high] or as a function handle with the form lfp_phases =
%   custom_phase_func(lfp_ts,lfp_sig) for custom phase estimation, for
%   example, by taking asymmetry into account
%   * slope_bnd: Default [].  Can be set to the slope bounds for anglereg.
%
%   RETURNS
%   * FIELD_INDEX: If method is 'grid', then the bins of the rate map are
%   percentile ranked between 0 and 1. Thus, a field index value of 0.5
%   indicates that 50% of the bins have a lower rate that the bin in
%   question. If method is 'place', then teh bins of the rate map are
%   linearly normalized between 0 and 1.
%
% This code has been freely distributed by the authors. If used or
% modified, we would appreciate it if you cited our paper:
% Climer, J. R., Newman, E. L. and Hasselmo, M. E. (2013), Phase coding by 
%   grid cells in unconstrained environments: two-dimensional phase 
%   precession. European Journal of Neuroscience, 38: 2526–2541. doi: 
%   10.1111/ejn.12256
%
% RELEASE NOTES
%   v1.0 2014-10-15 Release (Jason Climer, jason.r.climer@gmail.com)
%
% This file is part of pass_index.
%
% Copyright (c) 2014, Trustees of Boston University
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
%
% 1. Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in the 
% documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER 
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% Parse inputs
p = pass_index_parser(varargin{:});
for i = fields(p.Results)'
    eval([i{1} ' = p.Results.' i{1} ';']);
end
p.addParamValue('get_map',false);
p.parse(varargin{:});
get_map = p.Results.get_map;

% Make rate map & find positions of locations
[ map,~,occupancy,pos_loc,~ ] = rate_map( varargin{:} );
pos_loc(pos_loc==0) = NaN;

% Fill in unsampled holes
for d=1:size(pos,2)
    while any(isnan(pos_loc(:,d)))
        i = find(isnan(pos_loc(:,d)),1);
        if i==1
            j = find(~isnan(pos_loc(i:end,d)),1)+i-2;
            pos_loc(i:j,d) = pos_loc(j+1,d);
        elseif isempty(find(~isnan(pos_loc(i:end,d)),1))
            pos_loc(i:end,d) = pos_loc(i-1,d);
        else
            j = find(~isnan(pos_loc(i:end,d)),1)+i-2;
            temp = linspace(pos_loc(i-1,d),pos_loc(j+1,d),j-i+3);
            pos_loc(i:j,d) = temp(2:end-1);
        end
    end
end

% Convert locations to indicies
pos_loc = round(pos_loc);

% Addition to remove locations outside dim of map - Ryan H 2019-may-11
pos(size(map,1)<pos_loc(:,1) | size(map,2)<pos_loc(:,2),:)=[];
pos_loc(size(map,1)<pos_loc(:,1) | size(map,2)<pos_loc(:,2),:)=[];

temp = NaN(size(pos,1),1);
i = ~any(pos_loc==0,2);
pos_loc = cellfun(@(x)pos_loc(i,x),num2cell(1:size(pos,2)),'UniformOutput',false);
temp(i) = sub2ind(size(map),pos_loc{:});
pos_loc = temp;

% Make field index map
fi_map = map;
fi_map(occupancy==0) = NaN;

% Normalization methods
switch method
    case 'grid'
        % Normalize by percentile rate
        [fi_dist,i] = sort(fi_map(:));
        if any(isnan(fi_dist(:)))
            k = 1:find(isnan(fi_dist),1);
        else
            k = 1:numel(fi_dist);
        end
        fi_dist(k) = linspace(0,1,numel(k));
        fi_map(i) = fi_dist;
    case 'place'
        % Normalize to peak
        fi_map = (fi_map-nanmin(fi_map(:)))/range(fi_map(:));
end

if (get_map)
    field_index = fi_map;
else
    field_index = fi_map(pos_loc);
    field_index = interp1(pos_ts(~isnan(field_index)),field_index(~isnan(field_index)),pos_ts,'linear','extrap');
end