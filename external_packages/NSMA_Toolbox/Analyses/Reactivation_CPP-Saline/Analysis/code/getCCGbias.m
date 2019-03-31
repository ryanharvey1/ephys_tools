function [bias, pre, post] = getCCGBias(ccgca, norm, width)
%
%  bias = getCCGBias(ccgca, norm)
%
%  get bias estimate for all cell-pairs in a crosscorrelogram-cellarray ccgca.
%
% INPUT:
%   ccgca =   cross-correlogram cell array ccgca{1..npairs} of cross-correlograms (nbins+1 vectors);
%   norm  =   flag: 0 ... bias is unnormalized:  B = postcenter area - precenter area
%                   1 ... bias is normalized: B = (postcenter area - precenter area)/(postcenter area + precenter area)  
%   width =   number of bins in pre- and post-center areas (center bin is NOT included) 
% 
% OUTPUT:
%   bias = column vector of same length as input cell array ccgca of (un)normalized bias values; one for each cross-correlogram in ccgca.
%   pre = column vector of same length as input cell array ccgca of precenter area; one for each cross-correlogram in ccgca.
%   post = column vector of same length as input cell array ccgca of postcenter area; one for each cross-correlogram in ccgca.
%
% PL April 18, 2003

NPairs = length(ccgca);
pre = zeros(NPairs,1);
post = zeros(NPairs,1);
centerBin = 1 + floor(length(ccgca{1})/2);
if width > centerBin
    error('input argument width is greater than the half width of the cross-correlograms');
end

for ii = 1:NPairs
    pre(ii) =  sum(ccgca{ii}(centerBin-width-1:centerBin-1));
    post(ii) = sum(ccgca{ii}(centerBin+1:centerBin+width+1));        
end % for ii

if(norm)
    warning off;
    bias = (post - pre)./(post + pre);
    warning on;
else
    bias = post - pre;
end %if
            
