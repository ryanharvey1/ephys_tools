function hgamma = EMalgorithm(X,initial_hgamma)
%Input
%X:                 time series data
%initial_hgamma:    initial hyperparameters
%Output
%hgamma:            optimized hyperparameters
%
%version 0.1
%2010-3-31 Takeaki Shimokawa
%Copyright (c) 2010, Takeaki Shimokawa All rights reserved.

T = diff(X);    %ISIs
N = length(T);  %number of ISIs

hgamma = initial_hgamma;
for j = 1 : 100 %iteration of EM algorithm
    Y = KalmanFilter(X,hgamma);
    
    hgamma_lambda = 0;
    for i = 1 : N-1
        hgamma_lambda = hgamma_lambda + (...
            Y(2,i+1)+Y(2,i)-2.*Y(3,i)+(Y(1,i+1)-Y(1,i)).^2 ...
            )./T(i);
    end
    hgamma_lambda = sqrt(hgamma_lambda./(N-1));
    
    hgamma_kappa = 0;
    for i = 1 : N-1
        hgamma_kappa = hgamma_kappa + (...
            Y(5,i+1)+Y(5,i)-2.*Y(6,i)+(Y(4,i+1)-Y(4,i)).^2 ...
            )./T(i);
    end
    hgamma_kappa = sqrt(hgamma_kappa./(N-1));
    
    hgamma = [hgamma_lambda, hgamma_kappa];
end

