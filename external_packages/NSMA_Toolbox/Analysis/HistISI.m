function [H, binsUsed] = HistISI(TS, varargin)

% H = HistISI(TS)
%
% INPUTS:
%      TS = a single ts object
%
% OUTPUTS:
%      H = histogram of ISI
%      N = bin centers
%
% If no outputs are given, then plots the figure directly.
%
% ADR 1998
% version L5.1
% status PROMOTED

% v4.1 uses log10 instead of log
% v5.0 now uses msec not timestamps! 
% v5.1 changed xlabel

%--------------------
if ~isa(TS, 'ts'); error('Type Error: input is not a ts object.'); end
epsilon = 1e-100;

ISI = diff(Data(TS)/10) + epsilon;
maxLogISI = max(log10(ISI));
nBins = 200;
H = ndhist(log10(ISI)', nBins, 0, maxLogISI);
binsUsed = logspace(0,maxLogISI,nBins);

%-------------------
if nargout == 0
   plot(binsUsed, H); 
   xlabel('ISI, ms');
   set(gca, 'XScale', 'log');
   set(gca, 'YTick', max(H));
end

   
   