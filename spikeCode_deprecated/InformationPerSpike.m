function [IpS,IpSmat] = InformationPerSpike(rY,pX)
%
%  IpS = InformationPerSpike(rY,pX)
%
%  Mutual Information per spike after Gothard et al., J.Neuroscience, Jan.15,(1996) pp.823-835
%
% INPUT:
%     rY = cell array of Occ weighted rate histograms ( nBins x 1 column vectors)
%     pX = column vector of Occupancy PROBABILITIES (normalized to sum(pX) = 1)
% 
% OUTPUT:
%     IpS = 1 x nCells row vector of lower limit on mutual Information rate per Spike (bits/sec)
%

% R = cat(2,rY{:}); %removed this because my data is not arranged in a cell
% array, Ben C. 2014
R = rY;
[nBins,nCells] = size(R);

avR=pX'*R;
avRmat = kron(ones(nBins,1),avR); % copy row vector avR into matrix with nBins identical rows

relR=R./avRmat;
warning off;
log_relR = log2(relR);
warning on;
ij= find(isinf(log_relR));    % find -Inf's (log(0)) and replace with 0's
log_relR(ij) = 0;

pXmat = kron(pX,ones(1,nCells)); % copy column vector pX into matrix with nCells indentical columns
IpSmat = pXmat.*relR.*log_relR;
IpS = sum(IpSmat);

