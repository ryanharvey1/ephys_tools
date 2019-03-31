function [ev,ev_rev,rm_s1,rm_s2,rs1_s2] = ExplainedVariance(cs1,cm,cs2)
% computes the explained variance (partial correlation) from given vectors of correlations cs1, cm and cs2.
% [ev,ev_rev,rm_s1,rm_s2,rs1_s2] = ExplainedVariance(cs1,cm,cs2)
% 
% INPUT:
% the cs1,cm,cs2 are supposedly same-length vectors of corrlation coeffients (the so-called R-matrices)
% of each valid cell pair of a  sleep1,behavior and sleep2 subsession triplet; but they don't need to be.
% Any triplet of same length vectors holding random observations, can be subjected to the explained variance
% formula, which computes the variance of cs2 explained by cm when cs1 is taken into account.
%
% OUTPUT:
% ev       ... explained variance of cs2 by cm subtracting the correlations of (cs1,cs2)
% ev_rev   ... 'reverse' explained variance of cs1 by cm subtracting the correlations of (cs1,cs2)
%              (i.e. the roles of cs1 and cs2 are interchanged)
% rm_s1,rm_s2,rs1_s2 ... the r-coefficients of Maze-Sleep1, Maze-Sleep2 and Sleep1-Sleep2 entering into the partial correlation formula
%
% part of the Explained Variance Sleep analysis script suite
% Version 3.1
% PL Sep 2005
%
% Note: previous version 2.2 of this function computed a WRONG ev_rev!!!

% make sure inputs are column vectors
[nr,nc] = size(cs1);
if nr == 1, cs1 = cs1'; end
[nr,nc] = size(cm);
if nr == 1, cm = cm'; end
[nr,nc] = size(cs2);
if nr == 1, cs2 = cs2'; end

% need at least vectors of length 2
if  length(cs1) < 2
    ev = NaN;
    ev_rev = NaN;
    return;
end

% make sure input vectors have same length
if length(cs1) ~= length(cm) | length(cs1) ~= length(cs2)
    error('input vectors cs1, cm and cs2 must have same lenghts!');
end

% compute the correlations of the correlations
rm_s2 = diag(corrcoef(cs2,cm),1);
rm_s1 = diag(corrcoef(cs1,cm),1);
rs1_s2 = diag(corrcoef(cs1,cs2),1);

% compute the explained variance (EV) for m,s2|s1
norm = sqrt((1 - rm_s1^2) * (1 - rs1_s2^2));
ev =( (rm_s2 - rm_s1 * rs1_s2) / norm)^2;    % explained variance of m_s2|s1

% compute the reversed explained variance (EV) with role of S1 and S2 interchanged for m,s2|s1
norm_rev = sqrt((1 - rm_s2^2) * (1 - rs1_s2^2));
ev_rev =( (rm_s1 - rm_s2 * rs1_s2) / norm_rev)^2; % explained variance of m_s1|s2
