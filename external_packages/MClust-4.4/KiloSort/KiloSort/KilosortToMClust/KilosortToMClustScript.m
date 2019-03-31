% Example script to extract Spikes from KiloSort output and put into MClust
% .spikes file.

% Filename of the raw data file which was fed to Kilosort
rawFname = 'C:\Users\hasz\Documents\Data-temp\Promoted\R422\R422-2017-05-18\R422-2017-05-18-KilosortRaw-24TT.dat';

% Path containing the DIN .dat files (with the timestamp information)
timePath = 'C:\Users\hasz\Documents\Data-temp\PrePromoted\R422\R422-2017-05-18';

% Path containing Kilosort output .npy files
fpath = 'C:\Users\hasz\Documents\Data-temp\Promoted\R422\R422-2017-05-18';

% Filename of the channel map to use
chanMapFname = 'C:\Users\hasz\Documents\Data-temp\Promoted\R422\R422_channel_map_24TT.txt';

% Extract the spikes and save .spikes and .clu files
% (this could take 2-3 hrs)
ExtractSpikesFromKilosort(rawFname, timePath, fpath, chanMapFname);

% Now you can load the spikes files into MClust (and the .clu files)
% First add LoadIntanSpikes.m to the MClust\LoadingEngines directory
% And add feature_Peak15to25.m to the MClust\Features directory