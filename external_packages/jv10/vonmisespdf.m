function p = vonmisespdf (x, mu, K)
% VONMISESPDF Von Mises probability density function (pdf)
%   Y = VONMISESPDF(THETA,MU,K) returns the pdf of the Von Mises 
%   distribution with mean MU and concentration parameter K, 
%   evaluated at the values in THETA (given in radians).
%
%   --> www.paulbays.com

p = exp(K*cos(x-mu)) / (2*pi * besseli(0,K));

