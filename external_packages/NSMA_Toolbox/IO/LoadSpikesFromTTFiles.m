function S = LoadSpikesFromTTFiles(ttfiles)
%
%  S = LoadSpikesFromTTFiles(ttfiles)
%
%  returns a S Cell array of Spike ts objects (just like S = LoadSpikes(tfiles)) but reads through
% a list of ttfiles (.tt or TT*.dat ) rather than .t (spike) files 
%
% this is useful for analyses of uncut tetrode files as a whole
%
%  PL Feb 2008

for ii = 1:length(ttfiles)
    tsii = ReadTTt(ttfiles{ii});
    S{ii} = ts(tsii);
end

