function [B, LL] = JV10_fit (X, T, NT)
% JV10_FIT (X, T, NT)
%   Returns maximum likelihood parameters B for a mixture model describing 
%   recall responses X in terms of target T, non-target NT, and uniform 
%   responses. Inputs should be in radians, -PI <= X < PI. Fitting is based 
%   on an EM algorithm with multiple starting parameters.
%
%   B = JV10_FIT (X, T, NT) returns a vector [K pT pN pU], where K is 
%   the concentration parameter of a Von Mises distribution capturing 
%   response variability, pT is the probability of responding with the 
%   target value, pN the probability of responding with a non-target 
%   value, and pU the probability of responding "randomly". 
%
%   [B LL] = JV10_FIT (X, T, NT) additionally returns the log likelihood LL.
%
%   Ref: Bays PM, Catalao RFG & Husain M. The precision of visual working 
%   memory is set by allocation of a shared resource. Journal of Vision 
%   9(10): 7, 1-11 (2009) 
%
%   --> www.paulbays.com

if (nargin<2 || size(X,2)>1 || size(T,2)>1 || size(X,1)~=size(T,1) || nargin>2 && ~isempty(NT) && (size(NT,1)~=size(X,1) || size(NT,1)~=size(T,1)))     
    error('Input is not correctly dimensioned'); return; 
end

n = size(X,1); 

if (nargin<3) NT = zeros(n,0); nn = 0; else  nn = size(NT,2); end

% Starting parameters
K = [     1   10  100];
N = [  0.01  0.1  0.4];
U = [  0.01  0.1  0.4];

if nn==0, N = 0; end

LL = -inf; B = [NaN NaN NaN NaN];

warning('off','JV10_function:MaxIter');

% Parameter estimates

for i=1:length(K)
    for j=1:length(N)
        for k=1:length(U)
            [b ll] = JV10_function(X,T,NT,[K(i) 1-N(j)-U(k) N(j) U(k)]);
            if (ll>LL)
                LL = ll;
                B = b;
            end
        end
    end
end

warning('on','JV10_function:MaxIter');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright 2010 Paul Bays. This program is free software: you can     %
%   redistribute it and/or modify it under the terms of the GNU General  %
%   Public License as published by the Free Software Foundation.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%