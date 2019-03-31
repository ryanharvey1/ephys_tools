function [ DIC ] = computeDIC(histAng,hdTuning,OverallFR)
%Computes the directional information content on unsmoothed ratemap for
%directional data. Adapted from Langston et al 2010. by LB March 2018
%   
probOrient = histAng./sum(histAng); 
reIC = hdTuning./OverallFR;
log_IC = log2(reIC); log_IC(isinf(log_IC)) = 0;
ICi = probOrient.*reIC.*log_IC;
DIC = sum(ICi);
end

