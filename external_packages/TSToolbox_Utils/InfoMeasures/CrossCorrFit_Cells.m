function L = CrossCorrFit_Cells(qPre,qPost,sPost,binT,ep)

% [L,params] = CrossValidationCrossCorrPrediction_Cells(qPre,qPost,sPost,ep)
% 
% computes the loglikelihood prediction of an ensemble of target spike trains given another ensemble of (predictor) spike trains
% The prediction is built upon a Generalized Linear Model. The algorithm
% searches for the optimal gaussian window that best predict the spike
% trains of the target cells given another neuron.
% 
% INPUTS:
%     qPre:     a binned spike train matrix of the predictive cells (could be smoothed)
%     qPost:    the binned spike train of the target cells, same bins and same smoothing
%     sPost:    tsdArray of target cell spike trains
%     ep:       an intervalset objet of the epoch during which the prediction is built
%     
% OUTPUTS:
%     L:        information rate (bit/s)
%     params:   a 3D matrix giving, for each pair, the optimal delay and
%               s.d. of the integration window
%
%  Dependencies: TStoolbox, ComputeCrossCorrPrediction (and dependencies
%  therein)

% Copyright (C) 2016 Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


dQpre       = Data(qPre);
nbCpre      = size(dQpre,2);
dQpost      = Data(qPost);
nbCpost     = size(dQpost,2);

timeBins    = 2.^(0:10);
delayBins   = [-50:5:50];

L = NaN(nbCpost,nbCpre,length(timeBins),length(delayBins));
timeVec     = Range(qPre);
ratePost    = Rate(sPost,ep);      

warning off
for t=1:length(timeBins)
    for d=1:length(delayBins)
        bins        = (-5*timeBins(t)-abs(delayBins(d)):+5*timeBins(t)+abs(delayBins(d)));
        Y           = normpdf(bins,delayBins(d),timeBins(t));
        preTrain    = convn(dQpre,Y(:),'same');
        
        for preC=1:nbCpre
            preSm   = preTrain(:,preC); 
            
            for posC=1:nbCpost
                postTrain   = dQpost(:,posC);

                if any(dQpre(:,preC)>0) && any(dQpost(:,posC)>0)
                                                                            
                    b       = glmfit(zscore(preSm),postTrain,'poisson');
                    f       = exp(preSm*b(2)+b(1));
                    L(posC,preC,t,d)     = SpkTrainValuation(sPost{posC},tsd(timeVec,f/binT),ratePost(posC));
                end
            end                

        end
    end

end
warning on

fprintf('done\n')
L = L/tot_length(ep,'s');