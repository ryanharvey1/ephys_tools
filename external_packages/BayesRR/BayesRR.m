function [time, rate, regularity] = BayesRR(X)
%Program for estimating the instantaneous rate and regularity of occurrence
%
%Input
%X:             time series data
%Output
%time:
%rate:          firing rate
%regularity:    firing regularity
%
%Reference
%Estimating instantaneous irregularity of neuronal firing.
%T. Shimokawa and S. Shinomoto, Neural Computation (2009) 21:1931-1951. 
%
%version 0.1
%2010-3-31 Takeaki Shimokawa
%Copyright (c) 2010, Takeaki Shimokawa All rights reserved.

T = diff(X);    %ISIs
mu = mean(T);

hgamma = [1 .* mu.^-1.5, 0.1 .* mu.^-0.5];  %initial value of hypereparameters
hgamma = EMalgorithm(X,hgamma);
Y = KalmanFilter(X,hgamma);

A = X;
A(length(X)) = [];
B = X;
B(1) = [];

time = (A+B)/2;
rate = Y(1,:);
regularity = Y(4,:);

