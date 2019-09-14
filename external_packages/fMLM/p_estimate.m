function p = p_estimate(d, spikes, times, tol)
% spikes and times are matrices with dimensions: nd ny nx
% d is a vector of length nd
% estimate matrix p(i1, i2) = (sum_j spikes(j, i1, i2) )/ ( sum_j( d(j)*times(j, i1, i2) )
% returns [-1 -1; -1 -1] if there's a problem

[nd ny nx] = size(spikes);
if( nd ~= length(d) )
    fprintf(1, 'd and spikes do not have the right dimensions - d=%d spikes=%dx%dx%d \n',...
        length(d), size(spikes));
    p = [-1 -1; -1 -1];
else
    temp = repmat(d, [1 ny nx]);
    denom = squeeze(sum(temp.*times, 1));
    % if denom is 0 (most likely as times=0) we don't know what rate p to predict,
    % if denom < tol, make p small: p = spikes*tol
    denom( find(denom<tol) ) = 1/tol;
    p = squeeze(sum(spikes, 1))./denom;
end