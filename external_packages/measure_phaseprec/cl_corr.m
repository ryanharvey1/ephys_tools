function [circ_lin_corr pval slope_deg phi0_deg RR] = cl_corr(lin_x, circ_y_deg, min_slope, max_slope)
   
   % Function to (1) fit a line to circular-linear data and (2) determine
   % the circular-linear correlation coefficient
   %
   % Inputs : 
   %            lin_x - linear data
   %            circ_y_deg - circular data in degrees
   %            min_slope - minimum slope of fitted line
   %            max_slope - maximum slope of fitted line
   
   % Outputs:
   %            circ_lin_corr - circular-linear correlation coefficient
   %            pval - p-value of the correlation coefficient
   %            slope_deg - slope of the fitted line in deg
   %            phi0_deg - phase offset of the fitted line in deg
   
  circ_y_rad = (circ_y_deg ./ 180) * pi; % convert to rad
  [phi0, slope, RR] = cl_regression(lin_x, circ_y_rad, min_slope, max_slope); % fit line to data
  circ_x = mod(2 * pi * abs(slope) * lin_x, 2 * pi); % convert linear variable to circular one
  
  p_uniform = 0.5; % criterion for uniformity test
  
  % test for uniformity, and which measure then should be used:
  [pval_x] = circ_rtest(circ_x);
  [pval_y] = circ_rtest(circ_y_rad);
  if (pval_x > p_uniform) || (pval_y > p_uniform)
    [circ_lin_corr  pval] = circCorrUniform(circ_x', circ_y_rad');
  else
    [circ_lin_corr  pval] = circCorr(circ_x', circ_y_rad');
  end
  
  slope_deg = slope.*360; % slope [deg]
  phi0_deg = (phi0./pi).*180; % phase offset [deg]