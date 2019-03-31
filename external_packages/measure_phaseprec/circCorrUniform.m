function [rho p ts] = circCorrUniform(a,b)
%
%CIRCCORRUNIFORM computes a  circular correlation coefficient
%according to Jammmalamadaka 2001, page 177, equation 8.2.4, which  
%can deal with uniform distributions of a or b.   
%This function is equivalent to circCorrJammalamadaka2.m
%   
% Input
% a	angles (samples)
% b	angles (samples)
%
% Output
% rho   corr. coeff.
% p	significance probability
%
% Written by Richard Kempter Feb 13, 2008,
% to deal with uniform distributions of a or b (book on page 177 , ii)
% This function is an extension of the matlab function circCorr.m 
% by philipp berens
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens

if size(a,2) > size(a,1)
	a = a';
end

if size(b,2) > size(b,1)
	b = b';
end

n = length(a);
a_bar = circMean(a);
b_bar = circMean(b);

aplusb = a + b;
aminusb = a - b;

aplusb_bar =  circMean(aplusb);
aminusb_bar = circMean(aminusb);

R_aplusb =  circResLength(aplusb);
R_aminusb = circResLength(aminusb);


den = 2*sqrt(sum(sin(a - a_bar).^2) .* sum(sin(b - b_bar).^2));

rho =  n* (R_aminusb - R_aplusb) / den;	



%RK not sure wheter the next equations on the significance are still valid

l20 = mean(sin(a - a_bar).^2);
l02 = mean(sin(b - b_bar).^2);
l22 = mean((sin(a - a_bar).^2) .* (sin(b - b_bar).^2));

ts = sqrt((n * l20 * l02)/l22) * rho;
p = 2 * (1 - normcdf(abs(ts)));

