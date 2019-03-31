function [cv timeBound] = CV(rg,varargin);

% [cv timeBound] = CV(rg) computes the coefficient of variation of a spike
% train rg given in second. It computes std and mean of the difference between spike timing
% up to the 99% percentile of the cumulative histogram (given by
% timeBound) to prevent any divergence of the computed value for large
% ranges. 
% By default, any spike time difference exceeding 10 secdons will be
% discarded unless specified as a supplementary argument:
% [cv timeBound] = CV(rg,upperBound);
%
%Adrien Peyrache, 2013

if ~isempty(varargin)
    upperBound = varargin{1};
else
    upperBound = 10;
end

dt = diff(rg);
[h,b] = hist(dt,[0:0.002:upperBound]);
cs = cumsum(h)/sum(h);
ix = find(cs>0.95);
ix = ix(1);

timeBound = b(ix);

dt = dt(dt<timeBound);
cv = std(dt)/mean(dt);