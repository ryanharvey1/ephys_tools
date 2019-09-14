function field = find_2Dratemap(spikes, positions, boxcar)
% returns firing field with smoothing before dividing spikes by positions
% NB doesn't assume zero for unvisited positions: just ignores them.

b = ones(boxcar, boxcar);
c = ones(size(positions));
c( find(positions==0) ) = 0; 
denom = filter2(b, c);
denom(find(denom==0)) = NaN;

fpositions = filter2(b, positions);
fpositions = fpositions./denom;

fspikes = filter2(b, spikes);
fspikes = fspikes./denom;

field = fspikes./fpositions;

% set field = 0 in unoccupied locations
field(find(positions==0) ) = 0;