function WriteKiloSortChannelMap(chanMapIn, chanMapOut)
% Write a channel map file for Kilosort
%
% USAGE
%   WriteKiloSortChannelMap(chanMapIn, chanMapOut, fs);
%
% INPUT
%   chanMapIn - input channel map filename(in lab matlab format)
%   chanMapOut - filename to write KiloSort-formatted channel map file to
%
% OUTPUT
%   Writes the channel map file for KiloSort to chanMapOut
%
% June 2017
% Brendan Hasz
% haszx010@umn.edu
% David Redish Lab, University of Minnesota Twin Cities

% Load Intan -> TT/Si Probe channel map
eval(fileread(chanMapIn));

% Convert to KiloSort format
chanMap = 1:numel(channel_map);
connected = connected(:);
xcoords = xcoords(:)';
ycoords = ycoords(:)';
kcoords = shank(:)';

% Save KiloSort channel map as .mat file
save(chanMapOut,'chanMap','connected','xcoords','ycoords','kcoords','fs')

end