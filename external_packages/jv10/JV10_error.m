function [P B] = JV10_error(X, T)
% [P B] = JV10_ERROR (X, T)
%   Returns precision (P) and bias (B) measures for circular recall data. 
%   Inputs X and T are (nx1) vectors of responses and target values
%   respectively, in the range -PI <= X < PI. If T is not specified,
%   the default target value is 0.
%
%   Ref: Bays PM, Catalao RFG & Husain M. The precision of visual working 
%   memory is set by allocation of a shared resource. Journal of Vision 
%   9(10): 7, 1-11 (2009) 
%
%   --> www.paulbays.com

if nargin<2, T = zeros(size(X)); end    
if (any(abs(X)>pi) | any(abs(T)>pi)), error('Input values must be in radians, range -PI to PI'); return; end
if (any(size(X)~=size(T))), error('Inputs must have the same dimensions'); return; end
if size(X,1)==1, X = X'; T = T'; end

E = wrap(X-T); % error

% Precision
N = size(X,1);
x = logspace(-2,2,100); P0 = trapz(x,N./(sqrt(x).*exp(x+N*exp(-x)))); % Expected precision under uniform distribution
P = 1./cstd(E) - P0; % Corrected precision

% Bias
B = cmean(E);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright 2010 Paul Bays. This program is free software: you can     %
%   redistribute it and/or modify it under the terms of the GNU General  %
%   Public License as published by the Free Software Foundation.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%