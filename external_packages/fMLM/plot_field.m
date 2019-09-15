function field = plot_field(spikes, positions, boxcar)
% plots place field given dim1 x dim2 arrays of spikes and positions in (rows, columns) 
% i.e. (y down, x across) format
% field = plot_field(spikes, positions, boxcar);
% for no smoothing use boxcar = 1 or:
% imagesc(spikes./positions)
% colorbar

% smoothing before - NB don't assume zero for unvisited positions: just ignore.
field = find_2Dratemap(spikes, positions, boxcar);

% set unnoccupied positions to peak rate in sampled locations to be shown as white
% field(positions == 0) = 0;
% peak_rate = max(max(field));
% field(positions == 0) = (1+1/16)*peak_rate;
% 
% colormap([colormap(viridis(255)); 1 1 1]);
% imagesc(field);
% colorbar;



% find location of peak - y is down!
% ofset [ x; y] = center of 1st bin, scale = bin width (cm). Peak location is in cm.
% [y,x]=find(ffield == peak_rate);
% if >1 equal peaks, report average location in bins and in cm.
% peak_location = [ ofset(1)+(mean(x)-1)*scale; ofset(2)+(mean(y)-1)*scale ];