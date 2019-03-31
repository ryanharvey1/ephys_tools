function  L = crossCorrModel_Intensity(preTrain,postTrain,preTest,postTest,binT,params)

% Computes the cross-validated log-likelihood of a spike train given a (predictor) spike train in a
% gaussian temporal window of a given s.d. and mean.
%
%  USAGE
%
%    L = crossCorrModel_Likelihood(preTrain,postTrain,preTest,postTest,binT,params)
%
%    preTrain       pre-synaptic unsmoothed binned spike train from training set
%    postTrain      post-synaptic unsmoothed binned spike train from training set
%    preTest        pre-synaptic unsmoothed binned spike train from test set
%    postTest       post-synaptic unsmoothed binned spike train from test set
%    binT           bn size (in seconds)
%    params         vector of mean (1st element) and s.d. (2nd element) of
%                   the gaussian window

%    Dependencies: SpkTrainLogLikelihood

% Copyright (C) 2016 Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

Y       = normpdf((-10*params(2):10*params(2)),params(1),params(2));
preSm   = convn(preTrain,Y(:),'same');
b       = glmfit(preSm,postTrain,'poisson');
preSm   = convn(preTest,Y(:),'same');
f       = exp(preSm*b(2)+b(1));
L       = -SpkTrainLogLikelihood(postTest,f*binT);