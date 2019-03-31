function [ s,b ] = anglereg( x, theta, bnds)
% ANGLEREG - Linear-circular regression
%   Performes a linear-circular regression (Kempter et al, 2012)
%
%   [ S,B ] = ANGLEREG(X,THETA)
%   [ S,B ] = ANGLEREG(X,THETA,TOLERANCE)
%   [ S,B ] = ANGLEREG(X,THETA,TOLERANCE,MAXITER)
%
%   ARGUMENTS
%   * x: A vector of linear values
%   * theta: A vector, the same size as X, containing circular values (in
%   radians)
%   * bnds (optional): A 2 element vector with the lower and upper bounds
%   of the slope (B)
%
%   RETURNS
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
% This file is part of pass_index. Parts of this file may be
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
if ~exist('bnds','var')
    bnds = [];
end

% format inputs
theta = mod(theta,2*pi);
if size(x,1)>size(x,2)
    x = x';
end
if size(theta,1)>size(theta,2)
    theta = theta';
end

X = [ones(size(x(:))) x(:)];



if isempty(bnds)% No bounds
    % Initial guess from least squares and rolling phase
    phi = fminbnd(@(phi)sum((mod(theta(:)+phi,2*pi)-X*(inv(X'*X)*X'*mod(theta(:)+phi,2*pi))).^2),0,2*pi,...
    optimset('Display','off','TolX',1e-10));
    b = inv(X'*X)*X'*mod(theta(:)+phi,2*pi);
    b(1) = b(1)-phi;
    n = length(x);
    
    % Use above initial guess...
    [s1,r1] = fminsearch(...
        @(s)-sqrt((1/n*sum(cos(theta-2*pi*s*x)))^2+(1/n*sum(sin(theta-2*pi*s*x)))^2),...
        b(2));
    % And the perpendicular slope....
    [s2,r2] = fminsearch(...
        @(s)-sqrt((1/n*sum(cos(theta-2*pi*s*x)))^2+(1/n*sum(sin(theta-2*pi*s*x)))^2),...
        -1/b(2)...
        ,optimset('Display','off','TolX',1e-100,'TolFun',1e-100));
    % And take the better one
    if -r1>-r2
        s = s1;
    else
        s = s2;
    end
else
    % Use the bounds the user has given
    s = fminbnd(@(s)-sqrt((1/n*sum(cos(theta-2*pi*s*x)))^2+(1/n*sum(sin(theta-2*pi*s*x)))^2),bnds(1),bnds(2));
end

% Slope
b = atan2(sum(sin(theta-2*pi*s*x)),sum(cos(theta-2*pi*s*x)));
