function field = plot_dirs(spikes, positions, boxcar)
% plots polar field given dim1 arrays of spikes and positions as column vectors
% field = plot_dirs(spikes, positions, boxcar);
% boxcar smoothing (1 = no smoothing). Requires image processing toolbox

% smoothing - just ignore non-occupied places. 
% imfilter requires image processing toolbox, takes filter as 2nd argument (unlike filter2)
field = find_circ_ratemap(spikes, positions, boxcar);

% theta = ([1:length(field)]' - 0.5).*(2*pi/length(field));
% 
% polar([theta; theta(1)], [field; field(1)]);


