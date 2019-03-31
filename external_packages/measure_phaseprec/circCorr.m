function [rho pval ts] = circCorr(alpha,beta)
%
% computes circular correlation coefficient
%
% input:
%	alpha	sample of angles in radians
%	beta	sample of angles in radians
%
% output:
%	rho		correlation coefficient
%	pval	significance probability
%
% references:
%   Topics in circular statistics, S.R. Jammalamadaka et al., p. 176
%
% PHB 3/19/2006 2:02PM 
%
% copyright (c) 2006 philipp berens
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens
% distributed under GPL with no liability
% http://www.gnu.org/copyleft/gpl.html



if size(alpha,2) > size(alpha,1)
	alpha = alpha';
end
if size(beta,2) > size(beta,1)
	beta = beta';
end

n = length(alpha);
alpha_bar = circMean(alpha);
beta_bar = circMean(beta);

num = sum(sin(alpha - alpha_bar) .* sin(beta - beta_bar));
den = sqrt(sum(sin(alpha - alpha_bar).^2) .* sum(sin(beta - beta_bar).^2));

rho = num / den;	% correlation coefficient

l20 = mean(sin(alpha - alpha_bar).^2);
l02 = mean(sin(beta - beta_bar).^2);
l22 = mean((sin(alpha - alpha_bar).^2) .* (sin(beta - beta_bar).^2));

ts = sqrt((n * l20 * l02)/l22) * rho;
pval = 2 * (1 - normcdf(abs(ts)));

