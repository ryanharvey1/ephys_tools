function d = d_estimate(p, spikes, times, tol)
% spikes and times are matrices with dimensions: nd ny nx
% p is a matrix of size ny nx (y down)
% estimate vector d(j) = (sum_i1i2 spikes(j, i1, i2) )/ ( sum_i1i2( p(i1, i2)*times(j, i1, i2) )    
% returns [-1 -1] if there's a problem

[nd ny nx] = size(spikes);
if( sum([ny nx] == size(p))~= 2 )
    fprintf(1, 'p and spikes do not have the right dimensions - p=%dx%d spikes=%dx%dx%d \n',...
        size(p), size(spikes));
    d = [-1 -1];
else
    temp = repmat(p, [1 1 nd]);
    temp = permute(temp, [3 1 2]);
    denom = squeeze(sum(sum(temp.*times, 2), 3));
    % if denom is 0 (most likely as times=0) we don't know what rate p to predict,
    % if denom < tol, make p small: p = spikes*tol
    denom( denom<tol ) = 1/tol;
    d = squeeze(sum(sum(spikes, 2), 3))./denom;
end