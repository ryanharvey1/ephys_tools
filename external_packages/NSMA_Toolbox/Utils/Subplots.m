function [nL, ppL] = Subplots(n)

% Subplots  Chooses dimensions of nicely laid out grid for n subplots to be arranged in
%
% [nL, ppL] = Subplots(n)
%
% INPUTS:
%       n = number of plots to lay out
% OUTPUTS:
%       nL = number of lines
%       ppL = plots per line
%
% returns a pairing nL, ppL so that n plots are nicely laid out
% in a grid.  Then you can call subplot(nL, ppL, i) for each plot.
%
% ADR 1998, version v4.0, last modified '98 by ADR

% status PROMOTED


ppL = floor(sqrt(n));
nL = ceil((n)/ppL);



