function [phi0,slope,R] = cl_regression(x,phase,min_slope,max_slope)
%
%CL_REGRESSION determines the best linear fit to data on the
%              surface of a cylinder

%   Written by Richard Kempter (October-5, 2007)
%  
%   INPUTs:
%   x:      real-valued vector of `linear' instances 
%           (e.g. places, attenuations, frequencies, etc.)
%   phase:  vector of phases at instances x in rad (values need NOT be
%           restricted to the interval [0,2pi[ )
%   [min_slope,max_slope]:  interval of slopes in which the best_slope 
%            is determined.
%            In contrast to linear regression, we MUST restrict the range of
%            slopes to some `useful' range determined by prior knowledge.
%             ATTENTION ! because this is a possible source of errors, 
%            in particular for low numbers of data points.
%          
%   OUTPUTs:
%   R:      mean resultant lenght of the residual distribution 
%            is a measure of the goodness of the fit
%            (also called vector strength).
%           Small R indicates a bad fit (worst case R=0)   
%           Large R indicates a good fit (best case R=1)
%   slope:  slope (at the maximum of R) 
%                 within the interval [min_slope, max_slope]
%   phi0:   initial phase (or phase offset) of a cyclic regression line; 
%           values of phi0 are always restricted to the interval [0,2pi]   

  
% Check for valid inputs of the function
if nargin<1
   error('MATLAB:cl_regression:NotEnoughInputs', 'Not enough input arguments.');
end

if length(x)~=length(phase)
      error('MATLAB:cl_regression:XYmismatch', 'The lengths of x and phase must match.');
end

if length(x)<2
      error('MATLAB:cl_regression', 'The length of x is too small: length(x)<2.');
end

if ~isnumeric(min_slope) || ~isscalar(min_slope) 
   error('MATLAB:cl_regression','The ''min_slope'' parameter must be a scalar');
end

if ~isnumeric(max_slope) || ~isscalar(max_slope) 
   error('MATLAB:cl_regression','The ''max_slope'' parameter must be a scalar');
end

% add also a criterion min_slope < max_slope


% We determine the value of the best `slope' using the MATLAB function fminbnd.
% Please note that we have a minus sign in front of the 
% function `goodness' because we need the maximum (not the minimum)

slope = fminbnd(@(opt_slope) -goodness(x,phase,opt_slope),min_slope,max_slope);

% Given the best `slope' we can explicitly calculate the goodness of the
% fit ...
R = goodness(x,phase,slope);

% ... and the phase offset phi0:
[g_cos,g_sin] = summed_vector(x,phase,slope);
phi0 = atan2(g_sin,g_cos);
if phi0<0,   % restrict phases to the interval [0,2pi[
  phi0 = phi0 + 2 * pi;
end



%-----------------------------------------------------
% Here are subfunctions that are visible only within this function.
% These functions are necessary for the optimization procedure

function value = goodness(x,phase,slope)
% Determines the goodness of a fit of data in the cyclic domain
% with respect to some linear function of a specific slope

[g_cos,g_sin] = summed_vector(x,phase,slope); 
value = sqrt(g_cos^2 + g_sin^2)/length(x);


%-----------------------------------------------------
function [g_cos,g_sin] = summed_vector(x,phase,slope)
% Returns the two components of a vector,
% which is the summed difference between a linear function and the data points.

ph = phase - 2*pi*x*slope;  % phase difference
g_cos = sum(cos(ph));       % sum up phase differences
g_sin = sum(sin(ph));
  
% above added for perfomrance. original code follows HS2013
% g_cos=0; 
% g_sin=0;
%
% for i=1:length(x),
%      ph = phase(i) - 2* pi * x(i) * slope;  %phase difference
%      g_cos = g_cos + cos(ph);   % sum up phase differences
%      g_sin = g_sin + sin(ph);
% end 