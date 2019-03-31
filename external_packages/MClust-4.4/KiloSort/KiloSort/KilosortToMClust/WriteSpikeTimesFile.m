function WriteSpikeTimesFile(fname, T)
% Write .t file with the spike times in T.
%
% USAGE
%   WriteSpikeTimesFile(fname, T);
%
% INPUT
%   fname - filename for the .t file
%   T - Vector of spike times (in seconds)
%
% OUTPUT
%   Writes a .t file containing spike times.
%
% FILE FORMAT
%   A .t file is a binary file of uint32s, each the time of a spike in 10s 
%   of microseconds.
%
% August 2017
% Brendan Hasz
% haszx010@umn.edu
% David Redish Lab, University of Minnesota Twin Cities

fid = fopen(fname, 'wb', 'b');
fwrite(fid, uint32(T*10000), 'uint32', 0, 'b');
fclose(fid);

end