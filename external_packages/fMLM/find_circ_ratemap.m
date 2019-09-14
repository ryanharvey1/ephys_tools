function field = find_circ_ratemap(spikes, positions, boxcar);
% smoothing - just ignore non-occupied places. 
% imfilter requires image processing toolbox, takes filter as 2nd argument (unlike filter2)

b = ones(boxcar);
c = ones(size(positions));
c( find(positions==0) ) = 0; 
denom = imfilter(c, b, 'circular');
denom(find(denom==0)) = NaN;

fpositions = imfilter(positions, b, 'circular');
fpositions = fpositions./denom;
fspikes = imfilter(spikes, b, 'circular');  
fspikes = fspikes./denom; 

field = fspikes./fpositions;