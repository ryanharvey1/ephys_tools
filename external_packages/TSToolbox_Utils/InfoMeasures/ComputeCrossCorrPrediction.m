function params = ComputeCrossCorrPrediction(preTrain,postTrain,sPost,binT,varargin)

% Finds the optimal gaussian window (s.d. and time-lag) that integrates pre-synaptic spikes to
% predict post synaptic spike train.
%
%  USAGE
%
%    weights = ComputeCrossCorrPrediction(preTrain,postTrain,sPost,rg,binT,r,binT)
%
%    preTrain       a tsd of pre-synaptic unsmoothed binned spike train from training set
%    postTrain      a tsd of post-synaptic unsmoothed binned spike train from training set
%    sPost          post-synaptic spikeTrain
%    rg             time vector of the binned spike trains
%    binT           bin size
%
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'nSamples'            number of samples (default = 40)
%     'nFeatures'           number of features (default = 3)
%     'peakSampleIndex'     position of peak (default = 16)
%    =========================================================================
%
%    Dependencies: crossCorrModel_Likelihood, SpkTrainLogLikelihood,
%    TStoolbox

% Copyright (C) 2016 Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

if ~isempty(varargin)
    params0 = varargin{1};
else
    %weights0 = [0 1];
    params0 = [0 1];
end

%Average firing rate
rg          = Range(preTrain);
preTrain    = Data(preTrain);
postTrain   = Data(postTrain);
r           = mean(postTrain)/binT;

% Loklikelihood function to maximize
L = @(x)-SpkTrainValuation(sPost,tsd(rg,crossCorrModel_Intensity(preTrain,postTrain,x)/binT),r);

problem.objective       = L;
problem.options         = optimoptions('fmincon','display','off');
problem.solver          = 'fmincon';
problem.x0              = [0 1];
problem.lb              = [-0.200/binT 1];
problem.ub              = [0.200 0.500]/binT;
%problem.lb              = 1;
%problem.ub              = 0.500/binT;
problem.ObjectiveLimit  = 1e-10;


%weights = fminunc(problem);
warning off
params = fmincon(problem);
warning on
keyboard
% Y       = normpdf((-5*params:5*params),0,params);
% [h,b]   = xcorr(convn(preTrain,Y(:),'same'),postTrain,500);
% mins    = LocalMinima(-h, params, 0);
% 
% [~,ixB] = min(abs(b));
% [~,ix]  = min(abs(mins-ixB));
% mins    = b(mins(ix(1)));
% 
% L = @(x)-SpkTrainValuation(sPost,tsd(rg,crossCorrModel_Intensity(preTrain,postTrain,[-mins x])/binT),r);
% 
% problem.objective       = L;
% problem.x0              = 1;
% warning off
% pfinal = fmincon(problem);
% warning on
% fprintf('%f %f\n',params,pfinal)
% params = zeros(2,1);
% params(1) = mins;
% params(2) = pfinal;
