function fit = loglikelihood(p, d, spikes, times, tol)
% spikes and times are matrices with dimensions: nd ny nx
% p is a matrix of size ny nx (y down), d is a vector of length nd
% find log likelihood of spikes(j, i1, i2) under model: expected(j, i1, i2) = p(i1, i2).d(j).t(j, i1, i2) 
% assuming Poisson noise, i.e.: p(n spikes) = expected^n exp(-n) / n!
% fit = sum_ji1i2( spikes(j, i1, i2)*log(expected(j, i1, i2)) - expected(j, i1, i2) - log( spikes(j, i1, i2)! )
% NB log(expected) is replaced by log(tol) where expected < tol (e.g. 0.000001) to avoid log(0) divergence.
% uses gammaln(n+1) = log(n!) 

[nd ny nx] = size(spikes);
if( sum(size(spikes) == [length(d) size(p)]) ~= 3 )
    fprintf(1, 'p, d and spikes do not have the right dimensions - d=%d p=%dx%d spikes=%dx%dx%d \n',...
        length(d), size(p), size(spikes));
    fit = 1;
else    
    temp1 = repmat(p, [1 1 nd]);
    temp1 = permute(temp1, [3 1 2]);
    temp2 = repmat(d, [1 ny nx]);
    expected = temp1.*temp2.*times;
    %
    % log diverges for expected number of spikes = 0, if less than tol, replace with tol. 
    %
    expected2 = expected;
    expected2( find(expected < tol) ) = tol;
    fit = sum(sum(sum( spikes.*log(expected2) - expected - reshape(gammaln(spikes + 1), [nd ny nx]))));
end