function [L,weights] = CrossValidationSpikePrediction_Ensembles(qPre,qPost,sPost,binT,ep)

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

% Parameters:
nbEp        = 5; %Number of cross-validations
nbParams = 2;
epV = regIntervals(ep,nbEp);
fprintf('Launching Cross-validated Cross-Corr info\n')

dQpre       = Data(qPre);
nbCpre      = size(dQpre,2);
dQpost      = Data(qPost);
nbCpost     = size(dQpost,2);

rg = Range(qPre);
dt = median(diff(rg));

L = NaN(nbEp,nbCpost,nbCpre);

weights = zeros(nbCpost,nbCpre,nbEp,nbParams);

for ii=1:nbEp
    fprintf('.')
    
    ix = (1:nbEp);
    ix(ii)=[];
    testEp = epV{ii};
    trainingEp = epV{ix(1)};
    for jj=1:nbEp-2
        trainingEp = union(trainingEp,epV{ix(jj)});
        trainingEp = mergeCloseIntervals(trainingEp,1);
    end
    
    for posC=1:nbCpost
        qpost       = tsd(Range(qPost),dQpost(:,posC));
        postTrain   = Data(Restrict(qpost,trainingEp));
        postTest    = Data(Restrict(qpost,testEp));
        
        for preC=1:nbCpre
            qpre        = tsd(Range(qPost),dQpre(:,preC));
            preTrain    = Data(Restrict(qpre,trainingEp));
            preTest     = Data(Restrict(qpre,testEp));
            rgTest      = Range(Restrict(qpre,testEp));
            
            if any(dQpost(:,posC)>0)
            
                % Determines the optimal gaussian integration window
                warning off
                s   = Restrict(sPost{posC},testEp);
                r   = rate(s,testEp);
                b   = glmfit(preTrain,postTrain,'poisson');
                f   = exp(preTest*b(2)+b(1));
                
                L(ii,posC,preC)     = SpkTrainValuation(s,tsd(rgTest,f/binT),r);
                warning on
            end
        end
    end
end

fprintf('done\n')
L = squeeze(nansum(L));
L = L/tot_length(ep,'s');
