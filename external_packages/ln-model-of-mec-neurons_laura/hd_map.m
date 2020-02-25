function [hd_grid,dirVec,direction] = hd_map(direction,nbins)
% input: 
%   - head direction vector from data.frames 
%   - number of directional bins 
% output: 
%   - hd_grid: length of direction vector by nbins, which is binned head direction over time.  
hd_grid = zeros(length(direction),nbins);
dirVec = 2*pi/nbins/2:2*pi/nbins:2*pi-2*pi/nbins/2;

for i = 1:numel(direction)
    
    % figure out the hd index
    [~, idx] = min(abs(direction(i)-dirVec));
    hd_grid(i,idx) = 1; 
  
end

return