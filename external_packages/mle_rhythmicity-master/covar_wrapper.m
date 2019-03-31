function [ out ] = covar_wrapper( phat, cov, covars_list, asa )
%COVAR_WRAPPER Converts covarying parameter inputs
%   This wrapper function converts the fit represented in phat and the
%   covarying rhythmicity problem represented in cov, covars_list.
%
%   [out] = covar_wrapper(phat, cov, covars_list);
%       Returns the wrapped output with the standard parameters for
%       MLE_RHYTHMICITY, with s dropped if s not in cov
%       (noskip=true in mle_rhythmicity). See doc rhythmicity_covar for
%       more information
%
%   [out] = covar_wrapper(phat,cov,covars_list,asa)
%       If asa=true, then phat is formatted for the amplitude of the rhythm
%       (a) and not the relative amplitude (r). See doc mle_rhythmicity for
%       more information.
%
%   See also mle_rhythmicity, cif_generator, rhythmicity_covar
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

PARAMS = {'tau','b','c','f','s','r'};% List of parameters - used for easier coding
PARAMS = PARAMS(ismember(PARAMS,fields(cov)));
if ~all(ismember(fields(cov),PARAMS))
    PARAMS = fields(cov);
end

if ~exist('asa','var')
    asa = false;
end

if asa
    phat = [phat(numel(cov.r)+1:end) phat(1:numel(cov.r))/(1-phat(numel(cov.r)+numel(cov.tau)+1))];
end

vect = @(x)x(:);

out = cell(size(fields(cov)));
out(cellfun(@(x)isequal(cov.(x),0),PARAMS)) = ...
    num2cell(phat(arrayfun(@(i)sum(cellfun(@(s)numel(cov.(s)),...
    PARAMS(1:i))),find(cellfun(@(x)isequal(cov.(x),0),PARAMS)))));

out(~cellfun(@(x)isequal(cov.(x),0),PARAMS)) = ...
    ...
    arrayfun(@(i)...
    covars_list(:,cov.(PARAMS{i})+1)*...
    vect(phat(sum(cellfun(@(s)numel(cov.(s)),PARAMS(1:i-1)))+(1:numel(cov.(PARAMS{i})))))...
    ,find(~cellfun(@(x)isequal(cov.(x),0),PARAMS))...
    ,'UniformOutput',false...
    );

end

