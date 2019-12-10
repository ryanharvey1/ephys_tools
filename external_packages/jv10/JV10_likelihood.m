function [LL L] = JV10_likelihood (B, X, T, NT)
% --> www.paulbays.com

if (nargin<3 || size(X,2)>1 || size(T,2)>1 || size(X,1)~=size(T,1) || nargin>3 && ~isempty(NT) && (size(NT,1)~=size(X,1) || size(NT,1)~=size(T,1)))     
    error('Input is not correctly dimensioned');
    return; 
end

if (B(1)<0 || any(B(2:4)<0) || any(B(2:4)>1) || abs(sum(B(2:4))-1) > 10^-4 || (nargin<4 && B(3)~=0))
    warning('Invalid model parameters');
    LL = NaN; L = nans(size(X));
    return;
end

n = size(X,1); 

E  = X-T; E = mod(E + pi, 2*pi) - pi;

Wt = B(2) * vonmisespdf(E,0,B(1));
Wu = B(4) * ones(n,1)/(2*pi);

if (nargin<4)
    L = sum([Wt Wu],2);
else
    nn = size(NT,2);
    NE = repmat(X,1,nn)-NT; NE = mod(NE + pi, 2*pi) - pi;
    Wn = B(3)/nn * vonmisespdf(NE,0,B(1));
    L = sum([Wt Wn Wu],2);
end

LL = sum(log(L));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright 2010 Paul Bays. This program is free software: you can     %
%   redistribute it and/or modify it under the terms of the GNU General  %
%   Public License as published by the Free Software Foundation.         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%