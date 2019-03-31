function [ rho,p,s,b ] = kempter_lincirc( x,theta,varargin )
%KEMPTER_LINCIRC - Linear-circular correlation
%   Performes a linear-circular correlation (Kempter et al, 2012)
%
%   [ RHO,P,S,B ] = kempter_lincirc(X,THETA)
%   [ RHO,P,S,B ] = kempter_lincirc(X,THETA,S,B)
%
%   ARGUMENTS
%   * X: A vector of linear values
%   * THETA: A vector, the same size as X, containing circular values (in
%   radians)
%   * S: The presumed best slope in cycles/unit. If not entered, uses
%   anglereg
%   * B: The presumed best phase shift in radians. If not entered, uses
%   anglereg
%
%   RETURNS
%   * rho: The linear-circular correlation coefficient
%   * p: The statical significance of the correlation
%   * S: The slope in cycles per unit. To convert to radians, multiply by
%   2*pi
%   * B: The phase offset in radians.
%
% This code has been freely distributed by the authors. If used or
% modified, we would appreciate it if you cited our paper:
% Climer, J. R., Newman, E. L. and Hasselmo, M. E. (2013), Phase coding by 
%   grid cells in unconstrained environments: two-dimensional phase 
%   precession. European Journal of Neuroscience, 38: 2526–2541. doi: 
%   10.1111/ejn.12256
%
% RELEASE NOTES
%   v1.0 2014-10-15 Release (Jason Climer, jason.r.climer@gmail.com)
%
% This file is part of pass_index. All or part of this file may be
% considered derivative.
%
% Copyright (c) 2014, Trustees of Boston University
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
%
% 1. Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright 
% notice, this list of conditions and the following disclaimer in the 
% documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
% TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER 
% OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% Parse inputs
if nargin==1, s=varargin{1}; varargin = {}; end;
if strcmp(class(x),'single')
    x = double(x);
end

goods = ~isnan(x)&~isnan(theta);
x = x(goods);
theta = theta(goods);

% Angular regression
if ~exist('s','var')
    [s,b]=anglereg(x,theta,varargin{:});
end

n = length(x);

% Calculate rho
phi = mod(s*x,2*pi);
theta = mod(theta,2*pi);
phi_ = angle(sum(exp(1i*phi))/n);
theta_ = angle(sum(exp(1i*theta))/n);
rho = abs(sum(sin(theta-theta_).*sin(phi-phi_))/...
    sqrt(sum(sin(theta-theta_).^2)*sum(sin(phi-phi_).^2)))*sign(s);

% Calculate p
lambda = @(i,j)n^-1*sum((sin(phi-phi_).^i).*(sin(theta-theta_).^j));
z = rho*sqrt(n*lambda(2,0)*lambda(0,2)/lambda(2,2));
p = 1-erf(abs(z)/sqrt(2));

end

