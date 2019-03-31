function R = circResLength(alpha, w)
%
% computes mean resultant length for circular data
%
% input:
%	alpha	sample of angles in radians
%	[w		weightings in case of binned angle data]
%
% output:
%	R		mean resultant length
%
% PHB 3/19/2006 1:38PM
%
% references:
%   Statistical analysis of circular data, N.I. Fisher
%   Topics in circular statistics, S.R. Jammalamadaka et al. 
%
% copyright (c) 2006 philipp berens
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens
% distributed under GPL with no liability
% http://www.gnu.org/copyleft/gpl.html

if size(alpha,2) > size(alpha,1)
	alpha = alpha';
end

if nargin<2
	w = ones(size(alpha));
end

r = w'*exp(1i*alpha);
R = abs(r)/length(alpha);

