function [ lambda ] = rhythmicity_pdf(varargin)
%RHYTHMIC_FUN Parametric likelihood function for lags
%   Returns the likelihood of a lag occuring
%
% INPUTS:
%   phat: If rhythmic_fun is passed only one argument, the second is the
%       vector form of the parameters (below). phat may also be a cell
%       array. If passed no arguments, phat = [log10(0.3) 0.2 log10(0.5) 10
%       0.2 0.8]
%
%   tau (log10(sec)): Exponential falloff rate of whole distribution
%   b: Baseline CI
%   c (log10(sec)): Falloff rate of rhythmicity magnitude
%   f (Hz): Frequency of rhythmicity
%   s: Skipping (optional)
%   v: 1/2 rise time of the refractory rise (optional)
%   r: Rhythmicity
%
%   tau-r must be scalars or vectors of the same size as lags.
%
% PARAMETERS:
%   max_lag [0.6 sec]: The width of the analysis window
%   lags [0:0.01:0.6 sec]: Times to analyze the CI
%   islog (false): Whether to return the log-likelihood
%
% RETURNS:
%   lambda: The likelihood for the selected lags
%
% Copywrite (c) 2015, Trustees of Boston University
% All rights reserved.
%
% This file is part of mle_rhythmicity revision 2.0. The last committed
% version of the previous revision is the SHA starting with 93862ac...
%
% This code has been freely distributed by the authors under the BSD
% license (http://opensource.org/licenses/BSD2-Clause). If used or
% modified, we would appreciate if you cited our paper:
%
% Climer JR, DiTullio R, Newman EL, Hasselmo ME, Eden UT. (2014),
% Examination of rhythmicity of extracellularly recorded neurons in the
% entorhinal cortex. Hippocampus, 25:460-473. doi: 10.1002/hipo.22383.

% PROCESS INPUT
if (nargin==0||(any(cellfun(@ischar,varargin))&&find(cellfun(@ischar,varargin),1)==1)) % If no parameter arguments are passed
    phat =  [log10(0.3) 0.2 log10(0.5) 10 0.2 0.8]; % Set to default
    varargin = {};
end

if (nargin==1||(any(cellfun(@ischar,varargin))&&find(cellfun(@ischar,varargin),1)==2))% If using phat
    phat = varargin{1};
    varargin = varargin(2:end);
end

if exist('phat','var')
    if ~iscell(phat), phat=num2cell(phat); end
else
    phat = varargin(1:min([find(cellfun(@ischar,varargin),1)-1 numel(varargin)]));
    varargin = varargin(min([find(cellfun(@ischar,varargin),1)-1 numel(varargin)])+1:end);
end

ip = inputParser;

ip.addParamValue('max_lag',0.6);% Lag max_lag
ip.addParamValue('lags',0:0.01:0.6);% Lags to evaluate
ip.addParamValue('islog',false);

% Parse input and flush to workspace
ip.parse(varargin{:});
for j=fields(ip.Results)'
    eval([j{1} '= ip.Results.' j{1} ';']);
end

% Get each parameter
tau = phat{1};
b = phat{2};
c = phat{3};
f = phat{4};
if numel(phat)>5
    s = phat{5};
else
    s = 0;
end
r = phat{end};

s(s==0) = realmin;

[cif,cif_int] = cif_generator('full');
if islog
   lambda = log(cif(lags,tau,b,c,f,s,r))-log(cif_int(max_lag,tau,b,c,f,s,r)); 
else
   lambda = cif(lags,tau,b,c,f,s,r)/cif_int(max_lag,tau,b,c,f,s,r);
end

end


